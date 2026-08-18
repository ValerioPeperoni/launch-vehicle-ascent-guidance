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


# ---------------------------------------------------------------------------
# 7. Ciclo 10: correzione del target dello Stadio 2 (quota + verticalita'
#    invece di vx/vh con quota libera) e verifica di stabilita' orbitale.
# ---------------------------------------------------------------------------


def test_verifica_orbita_stabile_regressione_oracoli_noti():
    # Oracolo di regressione (STATUS.md, Ciclo 9 e Ciclo 10, addendum
    # "design finale consolidato", punto 3): valori gia' validati
    # indipendentemente due volte (dall'orchestratore, poi riconfermati
    # nell'indagine numerica pre-piano del Ciclo 10). Tolleranza stretta
    # (1 metro assoluto su grandezze dell'ordine di decine/centinaia di
    # km) -- un fallimento qui e' un bug nella formula, non un caso limite.

    # Stato vecchio target (Ciclo 7/9): quota libera, vx/vh vincolati.
    r_vecchio = val.verifica_orbita_stabile(148050.2, 7755.084, -8.241)
    assert r_vecchio["perigeo"] == pytest.approx(-62501.27, abs=1.0)  # ~= -62.50 km
    assert r_vecchio["apogeo"] == pytest.approx(148270.77, abs=1.0)  # ~= 148.27 km
    assert r_vecchio["perigeo_valido"] is False

    # Stato nuovo target (Ciclo 10): quota vincolata a 200 km esatti, vh=0.
    r_nuovo = val.verifica_orbita_stabile(200000.0, 7721.045, 0.0)
    assert r_nuovo["perigeo"] == pytest.approx(-22774.60, abs=1.0)  # ~= -22.77/-22.78 km
    assert r_nuovo["apogeo"] == pytest.approx(200000.0, abs=1.0)
    assert r_nuovo["perigeo_valido"] is False


def test_target_quota_verticale_raggiunto_entro_tolleranza_stretta(assemblaggio):
    # Il NUOVO obiettivo (Ciclo 10) vincola ATTIVAMENTE quota finale e
    # velocita' verticale finale (a differenza del vecchio target vx/vh,
    # che lasciava la quota libera). Verificato che il solutore converge
    # a precisione prossima a quella di macchina (residuo_h ~1e-6 m,
    # vh_finale ~1e-8 m/s, misurato in fase di sviluppo): tolleranza qui
    # fissata con ampio margine sopra il valore osservato (1 m, 1 mm/s),
    # non artificialmente stretta ne' larga.
    assert assemblaggio["h_finale"] == pytest.approx(assemblaggio["h_target"], abs=1.0)
    assert abs(assemblaggio["residuo_h"]) < 1.0  # metri
    assert abs(assemblaggio["vh_finale"]) < 1e-3  # m/s


def test_perigeo_non_valido_senza_bonus_valido_con_bonus_rotazione(assemblaggio):
    # Questi due controlli sono ESPLICITAMENTE collegati (non due casi
    # isolati): il secondo SPIEGA E CHIUDE il primo. Con lo stato finale
    # reale della traiettoria Falcon 9 di questo modulo, il modello
    # rigoroso del progetto (nessuna rotazione terrestre, CLAUDE.md,
    # riga "Dinamica") produce un perigeo leggermente sotto la superficie
    # (~-21/-23 km) -- NON un'orbita stabile in senso stretto. Sommando
    # il bonus di velocita' dovuto alla rotazione terrestre per un sito
    # equatoriale/subtropicale come Cape Canaveral (28.5°N,
    # BONUS_ROTAZIONE_TERRESTRE ~= 408.7 m/s, gia' quantificato nel
    # disclaimer di CLAUDE.md), lo stesso identico stato diventa
    # un'orbita ellittica stabile con perigeo sopra la superficie: il
    # deficit trovato senza bonus (~67 m/s in vx, vedi
    # ``residuo_vx``/CLAUDE.md) e' interamente coperto e superato dal
    # bonus (408.7 m/s > 67 m/s) -- non un limite aperto del modello, una
    # conferma quantitativa di coerenza (vedi STATUS.md, Ciclo 10).
    assert assemblaggio["orbita_senza_bonus"]["perigeo_valido"] is False
    assert assemblaggio["orbita_con_bonus"]["perigeo_valido"] is True

    # Il deficit di velocita' (senza bonus) e' interamente coperto dal
    # bonus di rotazione: verifica numerica diretta della relazione
    # causale appena descritta, non solo dei due booleani.
    deficit_vx = abs(assemblaggio["residuo_vx"])
    assert val.BONUS_ROTAZIONE_TERRESTRE > deficit_vx


def test_identita_chiusura_perdite_regge_con_la_nuova_traiettoria(perdite):
    # Verifica ESPLICITA (non assunta) che la tolleranza di chiusura di
    # ``scompone_perdite`` regga con la nuova forma di traiettoria
    # (Ciclo 10: A, B diversi, h(t) del segmento tangente lineare ora
    # converge monotonamente al target invece di salire e ridiscendere).
    # Misurato in fase di sviluppo: residuo_chiusura ~1.5e-4 m/s
    # (~1.9e-8 relativo) -- stesso ordine di grandezza del Ciclo 7/9, la
    # griglia di default di ``scompone_perdite`` (5000/20000/20000 punti)
    # NON necessita retaratura. Stessa tolleranza gia' usata da
    # ``test_identita_chiusura_perdite_tolleranza_stretta`` sopra,
    # ripetuta qui come check dedicato e commentato per il Ciclo 10.
    assert abs(perdite["residuo_chiusura"]) < 0.01
    assert abs(perdite["residuo_chiusura_relativo"]) < 1e-6


def test_risolvi_stadio_2_minimi_quadrati_ancora_coperta_direttamente(assemblaggio):
    # Ciclo 10: ``assembla_traiettoria`` non usa piu' di default
    # ``risolvi_stadio_2_minimi_quadrati`` (sostituita da
    # ``risolvi_stadio_2_target_quota``). La funzione vecchia resta pero'
    # nel modulo INVARIATA (vincolo esplicito del piano) e va quindi
    # continuare a esercitarla direttamente, non lasciarla come codice
    # morto non testato. Riusa lo stato di handoff Stadio1->Stadio2 REALE
    # (identico per costruzione: lo Stadio 1/gravity turn non e' toccato
    # da questo ciclo) e il vecchio target vx/vh, con la stessa
    # tolleranza attesa gia' documentata nel Ciclo 7 (~30-40 m/s di
    # scarto residuo, quota finale NON vincolata/libera).
    vx_target = val.velocita_orbitale_target()
    vh_target = 0.0
    soluzione = val.risolvi_stadio_2_minimi_quadrati(
        assemblaggio["m_ignizione_2_effettiva"],
        val.SPINTA_2,
        assemblaggio["mdot2"],
        assemblaggio["tf2"],
        assemblaggio["x1"],
        assemblaggio["h1"],
        assemblaggio["vx1"],
        assemblaggio["vh1"],
        vx_target,
        vh_target,
        val.A0_GUESS_STADIO2,
        val.B0_GUESS_STADIO2,
    )
    assert soluzione["coerenza_multistart"] is True
    assert soluzione["residuo_modulo"] == pytest.approx(34.4, abs=2.0)

    # La quota finale con QUESTO vecchio solutore resta libera (non
    # vincolata): con lo stato reale converge a ~148 km (ben lontana dal
    # target di 200 km), confermando esplicitamente perche' serviva la
    # correzione del Ciclo 10 (non solo dichiarato, verificato qui).
    h_finale_vecchio = soluzione["risultato"].y[1, -1]
    assert h_finale_vecchio == pytest.approx(148_050.2, rel=1e-3)
    assert h_finale_vecchio < 0.9 * val.H_TARGET
