"""Verifica della dinamica di ascesa verticale pura, singolo stadio (Step 2).

Parametri veicolo usati in questo modulo: vedi il commento
"# PROVVISORIO" sotto. Non sono dati di un lanciatore reale: quello e'
esplicitamente riservato allo Step 6 (CLAUDE.md, tabella vincoli, riga
"Caso di validazione"). Questi valori NON vanno confrontati con il
benchmark delta-v 9.1-10.0 km/s del CLAUDE.md (quel check e' per lo
Step 6, con stadio completo, guida e dati reali): qui, essendo
un'ascesa verticale pura senza gravity turn, ci si aspetta un delta-v
raggiunto nettamente inferiore al limite ideale di Tsiolkovsky (perdite
gravitazionali massime possibili, nessuna orbita raggiunta).
"""

import math

import numpy as np
import pytest

from lanciatore.costanti import CD, G0
from lanciatore.dinamica import derivate_stato, integra_ascesa_verticale

# PROVVISORIO — Step 2, sostituito da dati reali allo Step 6.
# Ordini di grandezza per un singolo stadio a propellente liquido
# kerolox, confrontabili con un primo stadio orbitale reale ma NON
# copiati da uno specifico lanciatore (vedi STATUS.md, Ciclo 2, punto
# 1.2 per la giustificazione completa di plausibilita' di ciascun
# valore).
M0 = 50_000.0  # kg, massa iniziale (con propellente)
M_VUOTO = 5_000.0  # kg, massa a vuoto -> frazione propellente 90%
SPINTA = 800_000.0  # N -> T/W0 = 800000/(50000*G0) ~= 1.63
ISP = 300.0  # s, tipico kerolox a livello del mare
DIAMETRO = 3.0  # m, ordine di grandezza primo stadio orbitale reale
AREA = math.pi * (DIAMETRO / 2) ** 2  # ~= 7.07 m^2

# Grandezze derivate, riportate qui per trasparenza (non hard-codate
# come costanti separate, vedi STATUS.md punto 1.2):
# mdot = SPINTA/(ISP*G0) ~= 271.9 kg/s
# t_burn = (M0 - M_VUOTO)/mdot ~= 165.5 s
# delta-v ideale Tsiolkovsky = ISP*G0*ln(M0/M_VUOTO) ~= 6774 m/s


@pytest.fixture(scope="module")
def risultato():
    """Integrazione eseguita una sola volta e riusata da tutti i test."""
    return integra_ascesa_verticale(M0, M_VUOTO, SPINTA, ISP, CD, AREA)


def test_condizioni_iniziali(risultato):
    # 6. h(0) == 0, v(0) == 0 esatti.
    h0, v0, m0 = risultato.y[:, 0]
    assert h0 == 0.0
    assert v0 == 0.0


def test_massa(risultato):
    # 1. m(0) == m0 esatto; m(t) strettamente decrescente durante il
    # bruciamento; valore finale == m_vuoto entro tolleranza numerica
    # coerente con l'atol usato dall'integratore (atol massa = 1e-3).
    m = risultato.y[2]
    assert m[0] == pytest.approx(M0)
    assert np.all(np.diff(m) < 0)
    assert m[-1] == pytest.approx(M_VUOTO, abs=1e-2)


def test_nessun_nan_o_inf(risultato):
    # 3. Nessun NaN/inf in nessuna componente della soluzione.
    for componente in risultato.y:
        assert np.all(np.isfinite(componente))


def test_monotonia_quota_e_velocita(risultato):
    # 2. h(t) strettamente crescente (conseguenza diretta di v > 0 dopo
    # il lift-off). v(t) strettamente crescente finche' la spinta netta
    # e' positiva, cioe' finche' T > D(t) + m(t)*g(h(t)). Il criterio e'
    # derivato direttamente dal segno di dv/dt calcolato con le stesse
    # equazioni/array che l'integratore produce (derivate_stato,
    # riusata da lanciatore/dinamica.py), NON da un ricalcolo
    # indipendente del drag qui nel test: cosi' un eventuale fallimento
    # e' inequivocabilmente un bug (nel drag o nei parametri), non un
    # caso limite da valutare caso per caso.
    h, v, m = risultato.y
    t = risultato.t

    assert np.all(np.diff(h) > 0)

    mdot = SPINTA / (ISP * G0)
    args = (SPINTA, ISP, CD, AREA, mdot, M_VUOTO)
    dv_dt = np.array(
        [derivate_stato(t[i], risultato.y[:, i], *args)[1] for i in range(len(t))]
    )

    # Con T/W0 ~= 1.63 la spinta netta deve restare positiva per
    # l'intera fase di bruciamento in questo caso di test. Se questo
    # assert fallisse sarebbe un bug da investigare (drag anomalo o
    # parametri non plausibili), non un caso da nascondere allargando
    # la tolleranza (regola di lavoro CLAUDE.md).
    assert np.all(dv_dt > 0)
    assert np.all(np.diff(v) > 0)


def test_confronto_tsiolkovsky(risultato):
    # 4. v_finale deve essere strettamente minore del delta-v ideale di
    # Tsiolkovsky (nessuna gravita', nessun drag): verifica solo la
    # direzione della disuguaglianza, non un valore target preciso.
    v_finale = risultato.y[1, -1]
    v_ideale = ISP * G0 * math.log(M0 / M_VUOTO)

    assert v_finale < v_ideale

    rapporto = v_finale / v_ideale
    # Ascesa verticale pura senza gravity turn: perdite gravitazionali
    # massime possibili (cos(0)=1 per l'intera durata) -> rapporto
    # atteso relativamente basso. Comportamento atteso di questo caso
    # di test degenere, non un difetto da correggere (vedi STATUS.md
    # punto 1.5.4). Se il test fallisse, questo commento e' dove
    # cercare il rapporto effettivo per confronto.
    assert 0.0 < rapporto < 1.0


def test_evento_di_terminazione(risultato):
    # 5. L'integrazione deve fermarsi per l'evento "fine propellente" E,
    # esplicitamente, l'evento "impatto suolo" non deve mai essersi
    # attivato (t_events vuoto per quell'evento) — non basta dedurlo
    # indirettamente dal fatto che l'altro evento e' popolato.
    t_events_fine_propellente, t_events_impatto_suolo = risultato.t_events

    assert risultato.status == 1  # un evento terminale ha fermato l'integrazione
    assert len(t_events_fine_propellente) == 1
    assert len(t_events_impatto_suolo) == 0
