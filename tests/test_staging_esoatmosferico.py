"""Verifica dello staging in fase esoatmosferica (Step 8).

Vedi STATUS.md, Ciclo 8 (piano + addendum reviewer) per il design
completo dello scenario. Nessun numero del piano e' hardcodato in questo
modulo: tutte le metriche (continuita', errore "senza ricalcolo",
convergenza dei coefficienti) sono ricalcolate qui in codice a partire
dal risultato di ``esegui_staging_esoatmosferico``, non copiate dal
report numerico del piano.
"""

import numpy as np
import pytest

from lanciatore.costanti import G0
from lanciatore.guida_esoatmosferica import risolvi_coefficienti_tangente_lineare
from lanciatore.staging_esoatmosferico import (
    TF2,
    esegui_staging_esoatmosferico,
)

# Tolleranza di residuo del root-finder, riusata (non allargata) da
# guida_esoatmosferica.risolvi_coefficienti_tangente_lineare, per il
# criterio 3 (A2,B2 raggiungono il target entro QUELLA tolleranza).
from lanciatore.guida_esoatmosferica import TOLLERANZA_RESIDUO_FSOLVE

# Soglia di errore relativo (criterio 2 del piano, "prova di
# invalidita'"): il piano riporta un errore relativo del 13.85% con la
# continuazione a tempo assoluto, "sopra la soglia 10%". Qui si verifica
# solo la direzione della disuguaglianza (>10%), il valore esatto e'
# ricalcolato in codice, non hardcodato.
SOGLIA_ERRORE_RELATIVO_SENZA_RICALCOLO = 0.10


@pytest.fixture(scope="module")
def scenario():
    """Esegue lo scenario una sola volta per modulo (integrazioni costose)."""
    return esegui_staging_esoatmosferico()


def test_continuita_stato_allo_staging(scenario):
    """x,h,vx,vh esatti tra fine Segmento1 e inizio Segmento2; massa -m_strut1."""
    stato_fine_seg1 = scenario["risultato_seg1"].y[:, -1]
    stato_inizio_seg2 = scenario["risultato_seg2"].y[:, 0]

    x_fine, h_fine, vx_fine, vh_fine, m_fine = stato_fine_seg1
    x_inizio, h_inizio, vx_inizio, vh_inizio, m_inizio = stato_inizio_seg2

    # Riporto diretto (nessuna conversione trigonometrica, a differenza
    # dello Step 5/7): x,h,vx,vh devono coincidere esattamente.
    assert x_inizio == x_fine
    assert h_inizio == h_fine
    assert vx_inizio == vx_fine
    assert vh_inizio == vh_fine

    # La massa differisce esattamente della massa strutturale espulsa.
    assert m_fine - m_inizio == pytest.approx(scenario["m_strut1"], abs=0.0)

    # Coerenza interna con i valori esposti dal dict di ritorno.
    assert x_fine == scenario["x_eff"]
    assert h_fine == scenario["h_eff"]
    assert vx_fine == scenario["vx_eff"]
    assert vh_fine == scenario["vh_eff"]
    assert m_fine == scenario["m_eff"]
    assert m_inizio == scenario["m_ignizione2"]


def test_evento_ode_genuino_non_cutoff_scriptato(scenario):
    """Il Segmento 1 si ferma per l'evento di massa, non per un t_span troncato a mano."""
    risultato_seg1 = scenario["risultato_seg1"]

    # status == 1 e' il codice solve_ivp per "un evento terminale ha
    # fermato l'integrazione" (diverso da 0 = t_span esaurito).
    assert risultato_seg1.status == 1

    # L'evento (unico, un solo evento passato a solve_ivp) deve avere
    # esattamente un istante registrato, e deve coincidere col tempo
    # finale della soluzione.
    assert len(risultato_seg1.t_events) == 1
    t_evento = risultato_seg1.t_events[0]
    assert len(t_evento) == 1
    assert t_evento[0] == pytest.approx(risultato_seg1.t[-1])

    # L'evento genuino si attiva PRIMA di tf1_nominale=40s (altrimenti
    # sarebbe indistinguibile da un semplice esaurimento di t_span):
    # verifica indipendente della massa raggiunta all'evento.
    assert risultato_seg1.t[-1] < 40.0

    # La massa finale deve coincidere (a tolleranza dell'integratore) con
    # la soglia fisica attesa: m0 - mdot1*durata_bruciamento_reale.
    m_soglia_attesa = scenario["risultato_seg1"].y[4, 0] - scenario["mdot1"] * (
        scenario["t_stage"]
    )
    assert risultato_seg1.y[4, -1] == pytest.approx(m_soglia_attesa, rel=1e-6)


def test_prova_invalidita_vecchi_coefficienti(scenario):
    """Errore euclideo normalizzato, continuazione a tempo assoluto, >10% del target."""
    risultato_senza_ricalcolo = scenario["risultato_senza_ricalcolo"]
    vx_finale = risultato_senza_ricalcolo.y[2, -1]
    vh_finale = risultato_senza_ricalcolo.y[3, -1]

    vx_target = scenario["vx_target"]
    vh_target = scenario["vh_target"]

    errore_euclideo = np.hypot(vx_finale - vx_target, vh_finale - vh_target)
    norma_target = np.hypot(vx_target, vh_target)
    errore_relativo = errore_euclideo / norma_target

    assert errore_relativo > SOGLIA_ERRORE_RELATIVO_SENZA_RICALCOLO

    # Sanity check sulla continuazione a tempo assoluto usata per questo
    # confronto: A0_continuato = A1 + B1*t_stage (non A1 invariato).
    assert scenario["A0_continuato"] == pytest.approx(
        scenario["A1"] + scenario["B1"] * scenario["t_stage"]
    )
    assert scenario["B0_continuato"] == scenario["B1"]


def test_coefficienti_ricalcolati_raggiungono_il_target(scenario):
    """A2,B2 (guess A0=A1,B0=B1) portano vx(tf2),vh(tf2) al target entro la tolleranza del risolutore."""
    risultato_seg2 = scenario["risultato_seg2"]
    vx_finale = risultato_seg2.y[2, -1]
    vh_finale = risultato_seg2.y[3, -1]

    residuo = np.array(
        [vx_finale - scenario["vx_target"], vh_finale - scenario["vh_target"]]
    )
    assert np.max(np.abs(residuo)) < TOLLERANZA_RESIDUO_FSOLVE * 10
    # Nota: la tolleranza di risolvi_coefficienti_tangente_lineare e'
    # verificata internamente sul RESIDUO DI FSOLVE (stesso sistema
    # risolto), qui si verifica il risultato dell'integrazione
    # indipendente via integra_tangente_lineare con un margine x10 per
    # assorbire l'errore di integrazione stesso (RK45 adattivo,
    # rtol/atol dell'integratore), senza allargare la tolleranza del
    # root-finder in se'.


def test_stabilita_regola_di_seeding(scenario):
    """Ripetere la risoluzione del Segmento 2 (stessa chiamata) da' A2,B2 identici."""
    ripetizioni = [
        risolvi_coefficienti_tangente_lineare(
            scenario["m_ignizione2"],
            scenario["spinta_seg2"],
            scenario["mdot2"],
            scenario["tf2"],
            scenario["vx_target"],
            scenario["vh_target"],
            scenario["A1"],
            scenario["B1"],
            x0=scenario["x_eff"],
            h0=scenario["h_eff"],
            vx0=scenario["vx_eff"],
            vh0=scenario["vh_eff"],
        )
        for _ in range(3)
    ]

    for A_rip, B_rip in ripetizioni:
        assert A_rip == pytest.approx(scenario["A2"])
        assert B_rip == pytest.approx(scenario["B2"])


def test_fenomeno_radici_multiple_guess_diverso(scenario):
    """Un guess diverso e non imparentato converge a una radice diversa ma valida."""
    A_diverso, B_diverso = risolvi_coefficienti_tangente_lineare(
        scenario["m_ignizione2"],
        scenario["spinta_seg2"],
        scenario["mdot2"],
        scenario["tf2"],
        scenario["vx_target"],
        scenario["vh_target"],
        0.5,
        -0.01,
        x0=scenario["x_eff"],
        h0=scenario["h_eff"],
        vx0=scenario["vx_eff"],
        vh0=scenario["vh_eff"],
    )

    # E' una radice DIVERSA da quella della regola di seeding A0=A1,B0=B1...
    assert not np.allclose(
        [A_diverso, B_diverso], [scenario["A2"], scenario["B2"]], atol=1e-3
    )

    # ...ma raggiunge comunque il target (verificato per integrazione,
    # non assunto): ricalcola la traiettoria con questi coefficienti e
    # confronta la velocita' finale col target.
    from lanciatore.guida_esoatmosferica import integra_tangente_lineare

    risultato_diverso = integra_tangente_lineare(
        scenario["m_ignizione2"],
        scenario["spinta_seg2"],
        scenario["mdot2"],
        A_diverso,
        B_diverso,
        scenario["tf2"],
        x0=scenario["x_eff"],
        h0=scenario["h_eff"],
        vx0=scenario["vx_eff"],
        vh0=scenario["vh_eff"],
    )
    vx_finale = risultato_diverso.y[2, -1]
    vh_finale = risultato_diverso.y[3, -1]
    residuo = np.hypot(
        vx_finale - scenario["vx_target"], vh_finale - scenario["vh_target"]
    )
    assert residuo < 1.0  # m/s, target raggiunto con ampio margine


def test_nessun_nan_inf_traiettoria_concatenata(scenario):
    """Nessun NaN/inf su Segmento1 + Segmento2 concatenati."""
    traiettoria_completa = np.concatenate(
        [scenario["risultato_seg1"].y, scenario["risultato_seg2"].y], axis=1
    )
    assert np.all(np.isfinite(traiettoria_completa))

    # Anche il risultato "senza ricalcolo" (usato solo per la prova di
    # invalidita', non fa parte della traiettoria reale del veicolo) deve
    # restare un numero finito, altrimenti il confronto del criterio 2
    # sarebbe privo di senso.
    assert np.all(np.isfinite(scenario["risultato_senza_ricalcolo"].y))


def test_tf2_coerente_con_tempo_residuo(scenario):
    """Sanity check: tf2 e' il tempo residuo nominale (40 - 15), non un valore scollegato."""
    assert TF2 == pytest.approx(40.0 - 15.0)
    assert scenario["tf2"] == TF2
