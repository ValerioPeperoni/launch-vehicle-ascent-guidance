"""Verifica del multistadio, fase atmosferica (Step 5).

Parametri di test marcati "# PROVVISORIO" (vedi STATUS.md, Ciclo 5,
addendum reviewer punto 5): tabella numerica a 2 stadi verificata
algebricamente dal reviewer, non dati di un lanciatore reale (riservato
allo Step 6).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo
modulo (riservato allo Step 6). I delta-v ideali per-stadio (Tsiolkovsky)
sono ricalcolati QUI con math.log a piena precisione, MAI hardcodati
come cifre scritte a mano (stesso principio gia' stabilito allo Step 4
per evitare errori di cancellazione catastrofica su un calcolo a mano,
vedi tests/test_guida_esoatmosferica.py).
"""

import math

import numpy as np
import pytest

from lanciatore.costanti import CD, G0
from lanciatore.staging import integra_multistadio_gravity_turn

# PROVVISORIO — Step 5, tabella numerica a 2 stadi CORRETTA (STATUS.md,
# Ciclo 5, punto 5 "Scoperta empirica del coder e correzione dei
# parametri di test"): m_prop_2/m_strut_2 ridotti da 8000/800 a 6000/600
# perche' con la tabella originale dell'addendum lo stadio 2 collassava
# per l'evento traiettoria_invalida (gamma->0) invece di completare il
# bruciamento nominale — vedi STATUS.md per la causa fisica (dgamma/dt
# non dipende dalla spinta, un secondo bruciamento troppo lungo sotto la
# stessa legge di guida fa collassare gamma). m0 = 51800 kg, payload
# implicito 1200 kg (invariato).
M0 = 51_800.0  # kg
STADIO_1 = {"m_prop": 40_000.0, "m_strut": 4_000.0, "spinta": 800_000.0, "isp": 300.0}
STADIO_2 = {"m_prop": 6_000.0, "m_strut": 600.0, "spinta": 150_000.0, "isp": 320.0}
STADI = [STADIO_1, STADIO_2]

# PROVVISORIO — Cd/area riusati dagli Step 2/3 (diametro 3.0 m).
DIAMETRO = 3.0  # m
AREA = math.pi * (DIAMETRO / 2) ** 2  # ~= 7.07 m^2

# PROVVISORIO — Step 3, stessi parametri di guida gia' usati allo Step 3.
V_KICK = 50.0  # m/s
KICK_ANGLE_DEG = 2.0  # gradi

# Payload atteso (implicito, mai un campo a se' nella struttura dati):
# m0 - m_prop_1 - m_strut_1 - m_prop_2 - m_strut_2, verificato a mano nel
# piano (STATUS.md, Ciclo 5, addendum reviewer punto 5).
PAYLOAD_ATTESO = 1_200.0  # kg


@pytest.fixture(scope="module")
def risultato_completo():
    """Integrazione multistadio completa, eseguita una sola volta e riusata."""
    return integra_multistadio_gravity_turn(
        M0, STADI, CD, AREA, V_KICK, KICK_ANGLE_DEG
    )


def test_1_continuita_stato_allo_staging(risultato_completo):
    # Criterio 1: x, h, v, gamma identici (entro la tolleranza numerica
    # dell'integratore) tra l'ultimo punto dello stadio 1 e il primo
    # punto dello stadio 2; la massa deve differire ESATTAMENTE (nessuna
    # tolleranza) di m_strut_1 (assegnazione diretta in staging.py, non
    # un risultato dell'integrazione).
    risultati = risultato_completo["risultati"]
    stato_fine_stadio_1 = risultati[0]["gravity_turn"].y[:, -1]
    stato_inizio_stadio_2 = risultati[1].y[:, 0]

    # x, h, v, gamma: indici 0..3, presi dallo stato effettivo in
    # entrambi i casi (nessun ricalcolo indipendente).
    assert stato_inizio_stadio_2[0] == pytest.approx(stato_fine_stadio_1[0])
    assert stato_inizio_stadio_2[1] == pytest.approx(stato_fine_stadio_1[1])
    assert stato_inizio_stadio_2[2] == pytest.approx(stato_fine_stadio_1[2])
    assert stato_inizio_stadio_2[3] == pytest.approx(stato_fine_stadio_1[3])

    # Massa: differenza ESATTA (non approssimata) di m_strut_1, perche'
    # e' un'assegnazione diretta in staging.py (massa_effettiva -
    # m_strut), non un valore prodotto dall'integratore.
    massa_fine_stadio_1 = stato_fine_stadio_1[4]
    massa_inizio_stadio_2 = stato_inizio_stadio_2[4]
    assert massa_fine_stadio_1 - massa_inizio_stadio_2 == STADIO_1["m_strut"]


def test_2_massa_effettiva_vs_nominale_al_burnout(risultato_completo):
    # Criterio (addendum 4): la massa EFFETTIVA all'evento di
    # fine-propellente di ogni stadio deve coincidere con la soglia
    # nominale m_vuoto_i entro l'atol dell'integratore (componente massa,
    # atol=1e-3 in staging.py/guida.py per lo stato a 5 componenti).
    riepilogo = risultato_completo["riepilogo"]

    for chiave in ("stadio_1", "stadio_2"):
        voce = riepilogo[chiave]
        assert voce["massa_burnout_effettiva"] == pytest.approx(
            voce["massa_burnout_nominale"], abs=1e-3
        )
        # Coerenza interna: differenza_massa deve essere esattamente
        # l'aritmetica dei due valori sopra.
        assert voce["differenza_massa"] == pytest.approx(
            voce["massa_burnout_effettiva"] - voce["massa_burnout_nominale"]
        )


def test_3_sanity_check_budget_massa_input():
    # Criterio 2/6 (rietichettato dall'addendum come sanity-check sui
    # dati di INPUT, vero per costruzione sui parametri di progetto, non
    # sulla simulazione): m0 - m_prop_1 - m_strut_1 - m_prop_2 -
    # m_strut_2 deve tornare esattamente al payload atteso.
    massa_residua = (
        M0
        - STADIO_1["m_prop"]
        - STADIO_1["m_strut"]
        - STADIO_2["m_prop"]
        - STADIO_2["m_strut"]
    )
    assert massa_residua == PAYLOAD_ATTESO


def test_4_confronto_tsiolkovsky_per_stadio(risultato_completo):
    # Criterio 3: delta-v ideale per-stadio (Tsiolkovsky, nessuna
    # gravita'/drag), ricalcolato qui con math.log (mai hardcodato).
    # m_ignizione_i e m_burnout_nominale_i sono le masse NOMINALI di
    # progetto (non gli effettivi della simulazione), coerenti con
    # l'identita' puramente algebrica del criterio.
    m_vuoto_1 = M0 - STADIO_1["m_prop"]
    m_ignizione_2 = m_vuoto_1 - STADIO_1["m_strut"]
    m_vuoto_2 = m_ignizione_2 - STADIO_2["m_prop"]

    dv_ideale_1 = STADIO_1["isp"] * G0 * math.log(M0 / m_vuoto_1)
    dv_ideale_2 = STADIO_2["isp"] * G0 * math.log(m_ignizione_2 / m_vuoto_2)
    dv_ideale_totale = dv_ideale_1 + dv_ideale_2

    risultati = risultato_completo["risultati"]
    v_finale_integrata = risultati[-1].y[2, -1]  # indice 2 = v in [x,h,v,gamma,m]

    assert v_finale_integrata < dv_ideale_totale


def test_5_nessun_evento_difensivo_attivato(risultato_completo):
    # Criterio 4: ogni evento di burnout deve essere "fine_propellente",
    # non un evento difensivo (impatto suolo, traiettoria invalida,
    # velocita' minima), in nessuno dei due stadi. Ordine degli eventi
    # in entrambi i segmenti: (fine_propellente, impatto_suolo,
    # traiettoria_invalida, velocita_minima) -> indici 1, 2, 3 devono
    # essere array vuoti.
    risultati = risultato_completo["risultati"]

    risultato_gt_1 = risultati[0]["gravity_turn"]
    risultato_2 = risultati[1]

    for risultato in (risultato_gt_1, risultato_2):
        for indice_evento_difensivo in (1, 2, 3):
            assert len(risultato.t_events[indice_evento_difensivo]) == 0


def test_6_nessun_nan_inf_su_traiettoria_concatenata(risultato_completo):
    # Criterio 5: nessun NaN/inf su tutta la traiettoria multistadio
    # concatenata (Fase A stadio 1 + Fase B stadio 1 + segmento stadio 2).
    risultati = risultato_completo["risultati"]

    fase_verticale_1 = risultati[0]["fase_verticale"]
    gravity_turn_1 = risultati[0]["gravity_turn"]
    stadio_2 = risultati[1]

    for risultato in (fase_verticale_1, gravity_turn_1, stadio_2):
        for componente in risultato.y:
            assert np.all(np.isfinite(componente))
