"""Verifica del caso di validazione con dati reali (Step 7, Ciclo 7, Fase B).

Esegue l'intera pipeline (``lanciatore.validazione.esegui_validazione``) UNA
SOLA VOLTA per modulo (fixture ``scope="module"``, la pipeline comporta
diverse integrazioni ODE, incluse le ri-integrazioni dense della
scomposizione delle perdite: ripeterla per ogni test singolarmente
sarebbe inutilmente lento) e verifica separatamente ciascun criterio
richiesto da STATUS.md, Ciclo 7, sezione "Criteri di verifica" del piano.
"""

import math

import numpy as np
import pytest

from lanciatore import validazione as val
from lanciatore.costanti import G0

# ---------------------------------------------------------------------------
# Fixture condivisa: l'intera pipeline gira una volta, tutti i test
# leggono lo stesso risultato (nessuna reintegrazione ripetuta).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def risultato_validazione():
    return val.esegui_validazione()


@pytest.fixture(scope="module")
def assemblaggio(risultato_validazione):
    return risultato_validazione["assemblaggio"]


@pytest.fixture(scope="module")
def perdite(risultato_validazione):
    return risultato_validazione["perdite"]


# ---------------------------------------------------------------------------
# 1. Dati veicolo citati (fonte nel docstring di modulo).
# ---------------------------------------------------------------------------


def test_fonte_dati_veicolo_citata_nel_docstring():
    # La fonte (Wikipedia, Falcon 9, consultata 2026-08-16, variante
    # espendibile) deve essere citata esplicitamente nel docstring di
    # modulo, non solo in STATUS.md -- stesso principio gia' applicato
    # in tutti gli step precedenti (es. Curtis Table 8.4 in atmosfera.py).
    docstring = val.__doc__
    assert "Wikipedia" in docstring
    assert "2026-08-16" in docstring
    assert "espendibile" in docstring
    assert "en.wikipedia.org/wiki/Falcon_9" in docstring


def test_valori_dati_veicolo_corretti():
    # Valori esatti richiesti dal piano (STATUS.md Ciclo 7, punto 2 e
    # addendum 3bis): verifica che nessuno sia stato alterato per errore.
    assert val.M_STRUT_1 == 25_600.0
    assert val.M_PROP_1 == 395_700.0
    assert val.SPINTA_1_VUOTO == 8_227_000.0
    assert val.ISP_1_VUOTO == 312.0
    assert val.ISP_1_SL == 283.0
    assert val.M_STRUT_2 == 3_900.0
    assert val.M_PROP_2 == 92_670.0
    assert val.SPINTA_2 == 981_000.0
    assert val.ISP_2 == 348.0
    assert val.PAYLOAD == 22_800.0
    assert val.DIAMETRO == 3.7
    assert val.AREA_RIFERIMENTO == pytest.approx(10.75, rel=1e-2)


def test_masse_nominali_bookkeeping():
    # STATUS.md Ciclo 7 punto 2: m0=540670 kg, m_vuoto_1=144970 kg,
    # m_ignizione_2=119370 kg, m_vuoto_2=26700 kg (verificati a mano
    # nell'addendum, qui ricontrollati contro il codice).
    assert val.M0_TOTALE == pytest.approx(540_670.0)
    assert val.M_VUOTO_NOMINALE_1 == pytest.approx(144_970.0)
    m_ignizione_2 = val.massa_ignizione_stadio_2(val.M_VUOTO_NOMINALE_1)
    assert m_ignizione_2 == pytest.approx(119_370.0)
    assert m_ignizione_2 - val.M_PROP_2 == pytest.approx(26_700.0)
    assert m_ignizione_2 - val.M_PROP_2 == pytest.approx(val.PAYLOAD + val.M_STRUT_2)


# ---------------------------------------------------------------------------
# 2. Continuita' Stadio1 -> Stadio2: (a) conversione v,gamma->vx,vh,
#    (b) sottrazione di m_strut_1 dalla massa EFFETTIVA, verificate
#    SEPARATAMENTE (addendum 3bis, punto 4).
# ---------------------------------------------------------------------------


def test_conversione_v_gamma_a_vx_vh_valori_noti():
    # v=100 m/s, gamma=30 gradi -> vx = 100*cos(30)=86.602540378...,
    # vh = 100*sin(30) = 50.0 esatto. Valori scelti per un confronto
    # analitico immediato (cos/sin di 30 gradi), indipendente
    # dall'assemblaggio completo.
    vx, vh = val.converti_v_gamma_a_vx_vh(100.0, math.radians(30.0))
    assert vx == pytest.approx(86.60254037844386, rel=1e-12)
    assert vh == pytest.approx(50.0, rel=1e-12)

    # Caso banale gamma=90 gradi (ascesa verticale pura): tutta la
    # velocita' e' verticale.
    vx2, vh2 = val.converti_v_gamma_a_vx_vh(200.0, math.radians(90.0))
    assert vx2 == pytest.approx(0.0, abs=1e-10)
    assert vh2 == pytest.approx(200.0, rel=1e-12)


def test_massa_ignizione_stadio_2_usa_massa_effettiva_non_nominale():
    # La funzione deve sottrarre m_strut1 da QUALUNQUE valore le venga
    # passato (la massa EFFETTIVA di fine Stadio 1), non essere agganciata
    # silenziosamente alla soglia nominale M_VUOTO_NOMINALE_1. Usiamo un
    # valore sintetico deliberatamente diverso dal nominale per dimostrarlo
    # in modo non tautologico.
    massa_effettiva_sintetica = 150_000.0  # != M_VUOTO_NOMINALE_1 (144970.0)
    risultato = val.massa_ignizione_stadio_2(massa_effettiva_sintetica)
    assert risultato == pytest.approx(150_000.0 - val.M_STRUT_1)
    assert risultato != pytest.approx(val.M_VUOTO_NOMINALE_1 - val.M_STRUT_1)


def test_continuita_stadio1_stadio2_nell_assemblaggio(assemblaggio):
    # Nell'assemblaggio reale, verifica che vx1/vh1 siano stati ottenuti
    # dalla conversione di v1/gamma1 (non ricalcolati altrove in modo
    # incoerente) e che la massa di ignizione Stadio 2 sia la massa
    # EFFETTIVA di fine Stadio 1 meno m_strut1.
    vx1_atteso, vh1_atteso = val.converti_v_gamma_a_vx_vh(assemblaggio["v1"], assemblaggio["gamma1"])
    assert assemblaggio["vx1"] == pytest.approx(vx1_atteso, rel=1e-12)
    assert assemblaggio["vh1"] == pytest.approx(vh1_atteso, rel=1e-12)

    m_ignizione_atteso = assemblaggio["m1_effettiva"] - val.M_STRUT_1
    assert assemblaggio["m_ignizione_2_effettiva"] == pytest.approx(m_ignizione_atteso, rel=1e-12)

    # Lo Stadio 1 deve essersi fermato per fine-propellente (non un evento
    # difensivo): se cosi' non fosse sarebbe un segnale da investigare
    # (stesso principio di staging.py, Step 5), non da ignorare.
    assert assemblaggio["stadio_1_fine_propellente"] is True


# ---------------------------------------------------------------------------
# 3. Delta-v ideale totale (vuoto) nel range 9.1-10.0 km/s, calcolato con
#    math.log in codice (CLAUDE.md, "Check di validazione obbligatorio").
# ---------------------------------------------------------------------------


def test_delta_v_ideale_totale_calcolato_con_math_log():
    # Ricalcolo indipendente con math.log direttamente nel test (non
    # importa la formula, la ricostruisce) per non essere tautologico.
    dv = val.calcola_delta_v_ideale()

    m0 = val.M0_TOTALE
    m_vuoto_1 = val.M_VUOTO_NOMINALE_1
    m_ignizione_2 = m_vuoto_1 - val.M_STRUT_1
    m_vuoto_2 = m_ignizione_2 - val.M_PROP_2

    dv1_atteso = val.ISP_1_VUOTO * G0 * math.log(m0 / m_vuoto_1)
    dv2_atteso = val.ISP_2 * G0 * math.log(m_ignizione_2 / m_vuoto_2)

    assert dv["dv1"] == pytest.approx(dv1_atteso, rel=1e-10)
    assert dv["dv2"] == pytest.approx(dv2_atteso, rel=1e-10)
    assert dv["dv_totale"] == pytest.approx(dv1_atteso + dv2_atteso, rel=1e-10)


def test_delta_v_ideale_totale_nel_range_claude_md():
    # CLAUDE.md, "Check di validazione obbligatorio": 9.1-10.0 km/s.
    dv = val.calcola_delta_v_ideale()
    assert 9100.0 <= dv["dv_totale"] <= 10_000.0


def test_delta_v_ideale_sensibilita_isp1_sl_riportato_esplicitamente():
    # Addendum 3bis punto 1: il numero con Isp1 SL deve essere disponibile
    # nel risultato (limite inferiore di sensibilita', non nascosto). Con
    # questi dati risulta ~8764 m/s, SOTTO il range -- e' il numero atteso
    # che motiva la scelta di usare Isp/spinta da vuoto per la traiettoria
    # principale (vedi docstring di modulo).
    dv_sl = val.calcola_delta_v_ideale(isp1=val.ISP_1_SL)
    assert dv_sl["dv_totale"] == pytest.approx(8763.8, abs=5.0)
    assert dv_sl["dv_totale"] < 9100.0  # sotto il range CLAUDE.md, come atteso e dichiarato


# ---------------------------------------------------------------------------
# 4. Identita' di chiusura della scomposizione delle perdite: i 4 termini
#    sommano a v_finale - v_iniziale entro tolleranza stretta (addendum
#    3bis punto 3, "check interno obbligatorio").
# ---------------------------------------------------------------------------


def test_identita_chiusura_perdite_tolleranza_stretta(perdite):
    # Verificato empiricamente (vedi docstring di lanciatore/validazione.py,
    # sezione "Nota tecnica importante sulla precisione del check di
    # chiusura") che con la ri-integrazione densa lo scarto di chiusura e'
    # ~1.5e-4 m/s (~2e-8 relativo). Tolleranza qui fissata con un margine
    # esplicito (~60x il valore osservato) sopra quel valore misurato, non
    # ne' artificialmente stretta ne' ampia per nascondere un problema:
    # 0.01 m/s assoluti su un v_finale di ~7755 m/s.
    assert abs(perdite["residuo_chiusura"]) < 0.01
    assert abs(perdite["residuo_chiusura_relativo"]) < 1e-6


def test_perdite_positive_e_finite(perdite):
    # Le tre perdite (gravita', drag, manovra) sono tutte >= 0 per
    # costruzione fisica (dv/dt viene sempre RIDOTTO da questi termini,
    # mai aumentato): un valore negativo significativo sarebbe un
    # sintomo di errore di segno, da investigare, non da accettare.
    assert perdite["perdita_gravita"] > 0.0
    assert perdite["perdita_drag"] > 0.0
    assert perdite["perdita_manovra"] >= 0.0  # puo' essere ~0 se theta=gamma quasi ovunque
    assert np.isfinite(perdite["perdita_gravita"])
    assert np.isfinite(perdite["perdita_drag"])
    assert np.isfinite(perdite["perdita_manovra"])


def test_budget_perdite_non_eccede_delta_v_ideale(perdite):
    # Vincolo fisico assoluto (indipendente dalla letteratura citata in
    # CLAUDE.md): la somma delle perdite non puo' MAI eccedere il delta-v
    # ideale disponibile, altrimenti v_finale sarebbe negativo, il che e'
    # fisicamente impossibile per una traiettoria che parte da v=0 e sale.
    somma_perdite = perdite["perdita_gravita"] + perdite["perdita_drag"] + perdite["perdita_manovra"]
    assert somma_perdite < perdite["dv_ideale_totale"]
    assert perdite["v_finale"] > 0.0


# ---------------------------------------------------------------------------
# 5. Scarto residuo dal target: piccolo in termini relativi (<1%) ma NON
#    forzato a zero -- si riporta cio' che emerge realmente (addendum
#    3bis punto 5).
# ---------------------------------------------------------------------------


def test_residuo_target_piccolo_ma_non_nullo(assemblaggio):
    residuo_relativo = assemblaggio["residuo_modulo"] / assemblaggio["vx_target"]
    # Piccolo in termini relativi (<1% del target), come atteso e
    # dichiarato esplicitamente nel piano (~0.4-0.5%, verificato).
    assert residuo_relativo < 0.01
    # NON deve essere (quasi) zero: se lo fosse, sarebbe il segno di una
    # convergenza artificialmente forzata, in contraddizione con la
    # decisione di design esplicita di NON usare il root-finder a
    # tolleranza stretta per questo problema (vedi docstring di modulo).
    assert assemblaggio["residuo_modulo"] > 1.0  # m/s, ben sopra il rumore numerico


def test_coerenza_multistart_ottimizzatore(assemblaggio):
    # Diagnostica di robustezza (non un gate bloccante nel codice, vedi
    # docstring di risolvi_stadio_2_minimi_quadrati): verificato che il
    # guess primario e quello perturbato convergono sullo stesso punto.
    # Se in futuro questo test fallisse, andrebbe investigato (nuovo
    # bacino di attrazione?), non semplicemente rimosso.
    assert assemblaggio["coerenza_multistart"] is True


# ---------------------------------------------------------------------------
# 6. Nessun NaN/inf in tutta la pipeline.
# ---------------------------------------------------------------------------


def test_nessun_nan_inf(assemblaggio, perdite):
    risultato_1 = assemblaggio["risultato_1"]
    for componente in risultato_1["fase_verticale"].y:
        assert np.all(np.isfinite(componente))
    for componente in risultato_1["gravity_turn"].y:
        assert np.all(np.isfinite(componente))
    for componente in assemblaggio["risultato_2"].y:
        assert np.all(np.isfinite(componente))

    for chiave in ("residuo_chiusura", "perdita_gravita", "perdita_drag", "perdita_manovra", "v_finale"):
        assert np.isfinite(perdite[chiave])

    for chiave in ("vx_finale", "vh_finale", "v_finale", "residuo_modulo", "A", "B"):
        assert np.isfinite(assemblaggio[chiave])
