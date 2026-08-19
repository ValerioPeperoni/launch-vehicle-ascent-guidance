"""Verifica del modulo di visualizzazione (Step 9, Ciclo 9).

Esegue la pipeline di validazione (``lanciatore.validazione.esegui_validazione``)
UNA SOLA VOLTA per modulo (fixture ``scope="module"``, stesso principio gia'
usato in ``tests/test_validazione.py``: la pipeline comporta diverse
integrazioni ODE, incluse le ri-integrazioni dense della scomposizione delle
perdite -- ripeterla per ogni test sarebbe inutilmente lento).

Criteri di verifica (STATUS.md, Ciclo 9, piano + addendum reviewer):
1. Lunghezza delle serie concatenate coerente coi dati sorgente.
2. Tempo monotono crescente attraverso le transizioni (continuita', non reset).
3. Marcatori di transizione coerenti con gli eventi reali.
4. Chiusura quantitativa (addendum punto 4): valori finali di v/h/m della
   serie concatenata contro ``assemblaggio['v_finale']``/valori finali di
   ``risultato_2``; massa ai due lati dello staging contro
   ``m1_effettiva``/``m_ignizione_2_effettiva``.
5. Le funzioni di grafico/animazione girano senza eccezioni e producono
   file non vuoti in ``output/``.
"""

import math
import os

import numpy as np
import pytest

from lanciatore import validazione as val
from lanciatore import visualizzazione as viz

# ---------------------------------------------------------------------------
# Fixture condivise.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def risultato_validazione():
    return val.esegui_validazione()


@pytest.fixture(scope="module")
def assemblaggio(risultato_validazione):
    return risultato_validazione["assemblaggio"]


@pytest.fixture(scope="module")
def serie(assemblaggio):
    return viz.estrai_serie_temporali(assemblaggio)


# ---------------------------------------------------------------------------
# 1. Lunghezza delle serie concatenate.
# ---------------------------------------------------------------------------


def test_lunghezza_serie_coerente_con_dati_sorgente(serie, assemblaggio):
    fa = assemblaggio["risultato_1"]["fase_verticale"]
    gt = assemblaggio["risultato_1"]["gravity_turn"]
    tl = assemblaggio["risultato_2"]

    lunghezza_attesa = len(fa.t) + len(gt.t) + len(tl.t)

    assert serie["n_fase_a"] == len(fa.t)
    assert serie["n_gravity_turn"] == len(gt.t)
    assert serie["n_tangente"] == len(tl.t)

    for chiave in ("t", "h", "v", "m", "x"):
        assert len(serie[chiave]) == lunghezza_attesa


# ---------------------------------------------------------------------------
# 2. Tempo monotono crescente attraverso le transizioni.
# ---------------------------------------------------------------------------


def test_tempo_monotono_crescente_attraverso_le_transizioni(serie):
    t = serie["t"]
    diff = np.diff(t)

    # Il tempo non deve MAI tornare indietro (nessun reset dell'orologio):
    # questo e' il criterio di continuita' richiesto dal piano.
    assert np.all(diff >= 0)

    # Le due transizioni (kick, staging) sono eventi ISTANTANEI per
    # costruzione (vedi validazione.py: la Fase A finisce e la Fase B
    # inizia allo stesso istante fisico, cosi' come il gravity turn finisce
    # e la tangente lineare inizia allo stesso istante fisico). Con
    # `offset_B = fase_verticale.t[-1]` e `gt.t[0] = 0`, il punto finale
    # della Fase A e il punto iniziale della Fase B cadono ESATTAMENTE
    # sullo stesso tempo assoluto (diff == 0 li', non un reset): stesso
    # discorso per la transizione staging -> tangente lineare. Ci si
    # aspettano quindi esattamente 2 diff nulli, esattamente alle 2
    # transizioni note -- non altrove (un diff nullo in un punto
    # imprevisto sarebbe un sintomo di bug, es. punti duplicati
    # dall'integratore).
    idx_diff_nulli = np.where(diff == 0.0)[0]
    idx_transizioni_attesi = [serie["n_fase_a"] - 1, serie["n_fase_a"] + serie["n_gravity_turn"] - 1]
    assert list(idx_diff_nulli) == idx_transizioni_attesi

    # Escludendo i 2 punti di transizione nota, il tempo e' strettamente
    # crescente ovunque.
    diff_senza_transizioni = np.delete(diff, idx_diff_nulli)
    assert np.all(diff_senza_transizioni > 0)


def test_offset_tempo_coerente_con_la_formula_dichiarata(serie, assemblaggio):
    fa = assemblaggio["risultato_1"]["fase_verticale"]
    gt = assemblaggio["risultato_1"]["gravity_turn"]

    offset_B_atteso = fa.t[-1]
    offset_tangente_atteso = offset_B_atteso + gt.t[-1]

    # Primo istante della Fase B (indice n_fase_a) = offset_B + gt.t[0] (=0).
    assert serie["t"][serie["n_fase_a"]] == pytest.approx(offset_B_atteso)
    # Primo istante della tangente lineare = offset_tangente + tl.t[0] (=0).
    idx_inizio_tangente = serie["n_fase_a"] + serie["n_gravity_turn"]
    assert serie["t"][idx_inizio_tangente] == pytest.approx(offset_tangente_atteso)


# ---------------------------------------------------------------------------
# 3. Marcatori di transizione coerenti con gli eventi reali.
# ---------------------------------------------------------------------------


def test_marcatori_transizione_coerenti_con_eventi_reali(serie, assemblaggio):
    fa = assemblaggio["risultato_1"]["fase_verticale"]
    gt = assemblaggio["risultato_1"]["gravity_turn"]

    idx_fine_fase_a, idx_fine_gravity_turn = serie["idx_transizioni"]

    # Indici attesi: ultimo punto di Fase A, ultimo punto di gravity turn
    # nell'array concatenato.
    assert idx_fine_fase_a == len(fa.t) - 1
    assert idx_fine_gravity_turn == len(fa.t) + len(gt.t) - 1

    # Il primo marcatore corrisponde alla quota/velocita' di kick (fine
    # Fase A verticale pura): la quota nel punto di transizione deve
    # coincidere con l'ultima quota della Fase A grezza.
    assert serie["h"][idx_fine_fase_a] == pytest.approx(fa.y[0, -1])
    assert serie["v"][idx_fine_fase_a] == pytest.approx(fa.y[1, -1])

    # Il secondo marcatore corrisponde allo staging (fine gravity turn):
    # la quota/velocita' nel punto di transizione deve coincidere con
    # l'ultimo punto grezzo del gravity turn (assemblaggio["h1"]/["v1"]).
    assert serie["h"][idx_fine_gravity_turn] == pytest.approx(assemblaggio["h1"])
    assert serie["v"][idx_fine_gravity_turn] == pytest.approx(assemblaggio["v1"])

    assert len(serie["t_transizioni"]) == 2
    assert len(serie["label_transizioni"]) == 2


# ---------------------------------------------------------------------------
# 4. Chiusura quantitativa (addendum reviewer, punto 4).
# ---------------------------------------------------------------------------


def test_chiusura_valori_finali_velocita_quota_massa(serie, assemblaggio):
    tl = assemblaggio["risultato_2"]

    h_finale_atteso = tl.y[1, -1]
    m_finale_atteso = tl.y[4, -1]
    v_finale_atteso = assemblaggio["v_finale"]

    assert serie["h"][-1] == pytest.approx(h_finale_atteso)
    assert serie["m"][-1] == pytest.approx(m_finale_atteso)
    assert serie["v"][-1] == pytest.approx(v_finale_atteso, rel=1e-9)


def test_chiusura_massa_ai_due_lati_dello_staging(serie, assemblaggio):
    # Nota di interpretazione (vedi report): l'addendum parla di "massa
    # INIZIALE di ciascun segmento" per entrambi i valori, ma
    # ``m1_effettiva`` e' per costruzione la massa FINALE del gravity turn
    # (fine bruciamento Stadio 1, vedi lanciatore/validazione.py,
    # `x1, h1, v1, gamma1, m1_effettiva = gt1.y[:, -1]`), non la massa
    # iniziale di un segmento. L'unica lettura coerente con i dati
    # disponibili e' verificare i due lati della discontinuita' di massa
    # allo staging: subito PRIMA (fine gravity turn) contro
    # m1_effettiva, e subito DOPO (inizio tangente lineare) contro
    # m_ignizione_2_effettiva -- che e' invece letteralmente la massa
    # iniziale del segmento tangente lineare.
    idx_fine_gravity_turn = serie["idx_transizioni"][1]
    idx_inizio_tangente = idx_fine_gravity_turn + 1

    assert serie["m"][idx_fine_gravity_turn] == pytest.approx(assemblaggio["m1_effettiva"])
    assert serie["m"][idx_inizio_tangente] == pytest.approx(assemblaggio["m_ignizione_2_effettiva"])

    # Il salto va PRESERVATO (nessuna interpolazione, addendum punto 5):
    # deve essere strettamente negativo (massa strutturale del primo
    # stadio espulsa) e di ordine di grandezza plausibile (~M_STRUT_1).
    salto_massa = serie["m"][idx_inizio_tangente] - serie["m"][idx_fine_gravity_turn]
    assert salto_massa < 0
    assert salto_massa == pytest.approx(-val.M_STRUT_1, rel=1e-6)


def test_massa_iniziale_segmento_tangente_coincide_con_m0_integrazione(serie, assemblaggio):
    # Controllo indipendente, letterale: la massa al primissimo istante
    # dell'integrazione dello Stadio 2 (`risultato_2.y[4, 0]`, il vero "m0"
    # passato a `integra_tangente_lineare`) deve coincidere esattamente
    # (stesso valore in memoria, nessuna operazione numerica in mezzo) con
    # `m_ignizione_2_effettiva`.
    tl = assemblaggio["risultato_2"]
    assert tl.y[4, 0] == pytest.approx(assemblaggio["m_ignizione_2_effettiva"])


def test_x_fase_a_sintetizzato_a_zero(serie):
    n_fase_a = serie["n_fase_a"]
    assert np.all(serie["x"][:n_fase_a] == 0.0)


# ---------------------------------------------------------------------------
# 5. Angolo di volo e serie interpolata (animazione avanzata): angolo
#    coerente con gamma/atan2 nei punti grezzi noti, staging preservato
#    come frame esatto del ricampionamento.
# ---------------------------------------------------------------------------


def test_angolo_volo_fase_a_costante_90_gradi(assemblaggio, serie):
    angolo = viz.angolo_volo_gradi(assemblaggio, serie)
    n_fase_a = serie["n_fase_a"]
    assert np.all(angolo[:n_fase_a] == 90.0)


def test_angolo_volo_gravity_turn_coincide_con_gamma_grezzo(assemblaggio, serie):
    gt = assemblaggio["risultato_1"]["gravity_turn"]
    angolo = viz.angolo_volo_gradi(assemblaggio, serie)
    n_fase_a = serie["n_fase_a"]
    n_gt = serie["n_gravity_turn"]
    assert np.allclose(angolo[n_fase_a : n_fase_a + n_gt], np.degrees(gt.y[3]))


def test_angolo_volo_tangente_lineare_coincide_con_atan2(assemblaggio, serie):
    tl = assemblaggio["risultato_2"]
    angolo = viz.angolo_volo_gradi(assemblaggio, serie)
    n_fase_a = serie["n_fase_a"]
    n_gt = serie["n_gravity_turn"]
    atteso = np.degrees(np.arctan2(tl.y[3], tl.y[2]))
    assert np.allclose(angolo[n_fase_a + n_gt :], atteso)


def test_angolo_volo_continuo_allo_staging(assemblaggio, serie):
    # La direzione della velocita' e' fisicamente continua allo staging
    # (solo la massa salta, vedi test_chiusura_massa_ai_due_lati_dello_staging):
    # l'angolo calcolato dal lato gravity turn (gamma) e dal lato tangente
    # lineare (atan2) devono quindi coincidere allo stesso istante.
    angolo = viz.angolo_volo_gradi(assemblaggio, serie)
    idx_fine_gt, idx_inizio_tl = serie["idx_transizioni"][1], serie["idx_transizioni"][1] + 1
    assert angolo[idx_fine_gt] == pytest.approx(angolo[idx_inizio_tl], abs=1e-6)


@pytest.fixture(scope="module")
def serie_interpolata(assemblaggio):
    return viz.estrai_serie_temporali_interpolata(assemblaggio, n_frame=220)


def test_serie_interpolata_istante_staging_esatto(serie_interpolata):
    i = serie_interpolata["idx_staging"]
    assert serie_interpolata["t"][i] == pytest.approx(serie_interpolata["t_staging"], abs=0.0)


def test_serie_interpolata_tempo_strettamente_crescente(serie_interpolata):
    # A differenza della serie grezza, qui l'istante di staging e' un
    # unico campione (non duplicato, vedi np.unique in
    # estrai_serie_temporali_interpolata): il tempo ricampionato deve
    # essere strettamente crescente ovunque, nessun diff nullo.
    assert np.all(np.diff(serie_interpolata["t"]) > 0)


def test_serie_interpolata_estremi_coincidono_con_serie_grezza(serie, serie_interpolata):
    assert serie_interpolata["t"][0] == pytest.approx(serie["t"][0])
    assert serie_interpolata["t"][-1] == pytest.approx(serie["t"][-1])
    assert serie_interpolata["h"][0] == pytest.approx(serie["h"][0])
    assert serie_interpolata["h"][-1] == pytest.approx(serie["h"][-1])
    assert serie_interpolata["v"][-1] == pytest.approx(serie["v"][-1])


def test_serie_interpolata_angolo_allo_staging_coincide_con_serie_grezza(assemblaggio, serie, serie_interpolata):
    idx_fine_gt = serie["idx_transizioni"][1]
    angolo_grezzo = viz.angolo_volo_gradi(assemblaggio, serie)
    i = serie_interpolata["idx_staging"]
    assert serie_interpolata["angolo_gradi"][i] == pytest.approx(angolo_grezzo[idx_fine_gt], abs=1e-6)


# ---------------------------------------------------------------------------
# 6. Grafici/animazione girano senza eccezioni e producono file non vuoti.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cartella_output_test(tmp_path_factory):
    return tmp_path_factory.mktemp("output_test")


def test_grafico_quota_velocita_tempo_produce_file_non_vuoto(serie, cartella_output_test):
    percorso = os.path.join(str(cartella_output_test), "quota_velocita_tempo.png")
    ritorno = viz.grafico_quota_velocita_tempo(serie, percorso)
    assert ritorno == percorso
    assert os.path.isfile(percorso)
    assert os.path.getsize(percorso) > 0


def test_grafico_traiettoria_produce_file_non_vuoto(serie, cartella_output_test):
    percorso = os.path.join(str(cartella_output_test), "traiettoria.png")
    ritorno = viz.grafico_traiettoria(serie, percorso)
    assert ritorno == percorso
    assert os.path.isfile(percorso)
    assert os.path.getsize(percorso) > 0


def test_anima_traiettoria_produce_file_non_vuoto(serie, cartella_output_test):
    percorso = os.path.join(str(cartella_output_test), "traiettoria.gif")
    ritorno = viz.anima_traiettoria(serie, percorso, n_frame=20, fps=10)
    assert ritorno == percorso
    assert os.path.isfile(percorso)
    assert os.path.getsize(percorso) > 0


def test_anima_traiettoria_avanzata_produce_file_non_vuoto(assemblaggio, cartella_output_test):
    # n_frame ridotto (non il default 220) solo per velocita' del test,
    # stesso principio gia' usato da test_anima_traiettoria_produce_file_non_vuoto
    # (n_frame=20) per l'animazione esistente.
    serie_interpolata_veloce = viz.estrai_serie_temporali_interpolata(assemblaggio, n_frame=20)
    percorso = os.path.join(str(cartella_output_test), "traiettoria_avanzata.gif")
    ritorno = viz.anima_traiettoria_avanzata(serie_interpolata_veloce, percorso, fps=10)
    assert ritorno == percorso
    assert os.path.isfile(percorso)
    assert os.path.getsize(percorso) > 0


def test_cartella_creata_automaticamente_se_non_esiste(serie, cartella_output_test):
    percorso = os.path.join(str(cartella_output_test), "sottocartella_nuova", "grafico.png")
    assert not os.path.isdir(os.path.dirname(percorso))
    viz.grafico_traiettoria(serie, percorso)
    assert os.path.isfile(percorso)


# ---------------------------------------------------------------------------
# Riepilogo testuale: riusa (non ricalcola) i numeri di validazione.py.
# ---------------------------------------------------------------------------


def test_riepilogo_testuale_contiene_i_numeri_di_validazione(risultato_validazione):
    testo = viz.riepilogo_testuale(risultato_validazione)
    assert isinstance(testo, str)

    dv_totale = risultato_validazione["delta_v_ideale"]["dv_totale"]
    v_finale = risultato_validazione["assemblaggio"]["v_finale"]

    # I valori numerici chiave devono comparire nel testo (stessa
    # precisione di formattazione usata da riepilogo_testuale, un
    # decimale), non ricalcolati qui in modo indipendente.
    assert f"{dv_totale:.1f}" in testo
    assert f"{v_finale:.1f}" in testo
