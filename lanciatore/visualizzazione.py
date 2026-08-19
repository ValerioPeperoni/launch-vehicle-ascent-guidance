"""Visualizzazione della traiettoria di ascesa (Step 9, Ciclo 9).

Converte i risultati grezzi di ``lanciatore.validazione.assembla_traiettoria``
(tre integrazioni ODE indipendenti, ciascuna con il proprio orologio interno
che riparte da t=0: Fase A ascesa verticale pura, Fase B gravity turn,
Stadio 2 tangente lineare) in serie temporali concatenate a TEMPO ASSOLUTO
CONTINUO, e produce grafici statici (quota/velocita'/massa nel tempo,
traiettoria nel piano (x,h)) e un'animazione GIF.

Nessuna nuova fisica in questo modulo: e' puro post-processing/plotting di
dati gia' calcolati e verificati allo Step 7. Riusa SENZA MODIFICHE
``lanciatore.validazione.esegui_validazione`` / ``assembla_traiettoria`` /
``calcola_delta_v_ideale`` / ``scompone_perdite`` (STATUS.md, Ciclo 9).

Formula esplicita dell'offset di tempo (STATUS.md Ciclo 9, addendum
reviewer, punto 1) -- i tre ``OdeResult`` grezzi sono indipendenti,
ciascuno riparte da t=0::

    offset_B = fase_verticale.t[-1]
    offset_tangente = offset_B + gravity_turn.t[-1]

Tabella di mapping stato -> grandezza, DIVERSA per ciascuna fase (STATUS.md
Ciclo 9, addendum reviewer, punto 2)::

    Fase A (verticale)    [h, v, m]             h=idx0  v=idx1 (diretto)             m=idx2  x=non esiste (sintetizzato a zero)
    Fase B (gravity turn) [x, h, v, gamma, m]   h=idx1  v=idx2 (diretto)             m=idx4  x=idx0
    Tangente lineare      [x, h, vx, vh, m]     h=idx1  v=sqrt(vx**2 + vh**2)        m=idx4  x=idx0

``x`` per la Fase A non esiste nei dati grezzi (ascesa verticale pura):
sintetizzato esplicitamente come ``np.zeros_like(fase_verticale.t)``
(addendum reviewer, punto 3) -- NON un IndexError/riuso implicito di
un'altra colonna.

Discontinuita' di massa allo staging (Stadio1 -> Stadio2, addendum
reviewer, punto 5): comportamento ATTESO e PRESERVATO deliberatamente
(nessuna interpolazione) -- concatenando gli array grezzi cosi' come sono,
la serie m(t) mostra correttamente un salto verso il basso nel punto di
transizione (massa strutturale del primo stadio espulsa allo staging).
"""

import os

import matplotlib

matplotlib.use("Agg")  # backend non interattivo: sicuro per generazione headless/test, nessuna finestra aperta
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.markers import MarkerStyle
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

from lanciatore import validazione as val

# Cartella di output di progetto (root del repo, non dentro lanciatore/):
# stesso trattamento di venv/ in .gitignore, artefatti rigenerabili.
_CARTELLA_OUTPUT_PROGETTO = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output"
)


def _assicura_cartella(percorso_file):
    """Crea la cartella del file di output se non esiste ancora."""
    cartella = os.path.dirname(os.path.abspath(percorso_file))
    if cartella and not os.path.isdir(cartella):
        os.makedirs(cartella, exist_ok=True)


def estrai_serie_temporali(assemblaggio):
    """Estrae e concatena le serie temporali delle 3 fasi a tempo assoluto continuo.

    Vedi docstring di modulo per la formula dell'offset di tempo e la
    tabella di mapping stato -> grandezza (diversa per ciascuna fase).

    Parametri
    ---------
    assemblaggio : dict
        Risultato di ``lanciatore.validazione.assembla_traiettoria`` (o del
        campo ``"assemblaggio"`` di ``esegui_validazione``).

    Ritorna
    -------
    dict
        ``{"t", "h", "v", "m", "x"}``: array numpy concatenati a tempo
        assoluto continuo. ``"idx_transizioni"``: indici (nell'array
        concatenato) dell'ULTIMO punto di ciascuna fase prima della
        transizione successiva (kick, staging). ``"t_transizioni"``:
        istanti assoluti corrispondenti. ``"label_transizioni"``:
        etichette leggibili. ``"n_fase_a"``, ``"n_gravity_turn"``,
        ``"n_tangente"``: lunghezze dei 3 segmenti sorgente (usate dai
        test per verificare la chiusura di lunghezza senza duplicare qui
        la logica di concatenazione).
    """
    fa = assemblaggio["risultato_1"]["fase_verticale"]  # OdeResult, stato [h, v, m]
    gt = assemblaggio["risultato_1"]["gravity_turn"]  # OdeResult, stato [x, h, v, gamma, m]
    tl = assemblaggio["risultato_2"]  # OdeResult, stato [x, h, vx, vh, m]

    offset_B = fa.t[-1]
    offset_tangente = offset_B + gt.t[-1]

    t_fa = fa.t
    t_gt = gt.t + offset_B
    t_tl = tl.t + offset_tangente

    # --- Fase A (ascesa verticale pura): [h, v, m], x sintetizzato a zero. ---
    h_fa, v_fa, m_fa = fa.y[0], fa.y[1], fa.y[2]
    x_fa = np.zeros_like(fa.t)

    # --- Fase B (gravity turn): [x, h, v, gamma, m], v gia' in modulo. ---
    x_gt, h_gt, v_gt, m_gt = gt.y[0], gt.y[1], gt.y[2], gt.y[4]

    # --- Tangente lineare: [x, h, vx, vh, m], v ricostruito dal modulo. ---
    x_tl, h_tl, vx_tl, vh_tl, m_tl = tl.y[0], tl.y[1], tl.y[2], tl.y[3], tl.y[4]
    v_tl = np.sqrt(vx_tl**2 + vh_tl**2)

    t = np.concatenate([t_fa, t_gt, t_tl])
    h = np.concatenate([h_fa, h_gt, h_tl])
    v = np.concatenate([v_fa, v_gt, v_tl])
    m = np.concatenate([m_fa, m_gt, m_tl])
    x = np.concatenate([x_fa, x_gt, x_tl])

    n_fase_a = len(t_fa)
    n_gravity_turn = len(t_gt)
    n_tangente = len(t_tl)

    # Ultimo indice di ciascuna delle prime 2 fasi = punto di transizione
    # (kick = fine Fase A, staging = fine gravity turn). Il salto di massa
    # allo staging e' PRESERVATO per costruzione (nessuna interpolazione,
    # vedi docstring di modulo): idx_fine_gravity_turn e idx_fine_gravity_turn+1
    # sono due punti temporalmente adiacenti (stesso istante assoluto, a
    # meno della precisione degli orologi separati) con massa diversa.
    idx_fine_fase_a = n_fase_a - 1
    idx_fine_gravity_turn = n_fase_a + n_gravity_turn - 1

    return {
        "t": t,
        "h": h,
        "v": v,
        "m": m,
        "x": x,
        "idx_transizioni": [idx_fine_fase_a, idx_fine_gravity_turn],
        "t_transizioni": [float(t[idx_fine_fase_a]), float(t[idx_fine_gravity_turn])],
        "label_transizioni": ["kick (fine Fase A)", "staging (fine gravity turn)"],
        "n_fase_a": n_fase_a,
        "n_gravity_turn": n_gravity_turn,
        "n_tangente": n_tangente,
    }


def angolo_volo_gradi(assemblaggio, serie):
    """Angolo di volo (rispetto all'orizzontale locale) per ciascun punto grezzo concatenato.

    Funzione isolata e testabile separatamente (stesso principio gia'
    seguito per ``lanciatore.validazione.converti_v_gamma_a_vx_vh``): non e'
    nuova fisica, e' una derivazione geometrica pura da stato gia'
    integrato, con la STESSA convenzione di ``gamma`` in ``guida.py``
    (90 gradi = verticale, gradi decrescenti verso l'orizzontale):

    - Fase A (ascesa verticale pura): ``vx=0`` per costruzione (nessuna
      componente orizzontale nel modello di questa fase) -> 90 gradi
      costanti.
    - Fase B (gravity turn): ``gamma`` e' gia' presente nello stato
      grezzo (``gt.y[3]``, radianti, stessa convenzione: vedi
      ``guida.py``, ``gamma0 = radians(90 - kick_angle_deg)``) -> solo
      conversione a gradi.
    - Tangente lineare: ricostruito con ``atan2(vh, vx)`` dalle
      componenti gia' presenti nello stato grezzo (``tl.y[2]``,
      ``tl.y[3]``) -- stessa convenzione (0=orizzontale, 90=verticale)
      per costruzione geometrica di ``atan2``.

    Parametri
    ---------
    assemblaggio : dict
        Stesso input di ``estrai_serie_temporali``.
    serie : dict
        Risultato di ``estrai_serie_temporali(assemblaggio)`` (riusato
        per ``n_fase_a``, non ricalcolato qui).

    Ritorna
    -------
    numpy.ndarray
        Angolo in gradi, stessa lunghezza/allineamento di ``serie["t"]``.
    """
    gt = assemblaggio["risultato_1"]["gravity_turn"]  # [x, h, v, gamma, m]
    tl = assemblaggio["risultato_2"]  # [x, h, vx, vh, m]

    angolo_fa = np.full(serie["n_fase_a"], 90.0)
    angolo_gt = np.degrees(gt.y[3])
    angolo_tl = np.degrees(np.arctan2(tl.y[3], tl.y[2]))

    return np.concatenate([angolo_fa, angolo_gt, angolo_tl])


def estrai_serie_temporali_interpolata(assemblaggio, n_frame=220):
    """Ricampiona la serie temporale a passo di tempo uniforme, per un'animazione fluida.

    Wrapper puramente di rendering attorno a ``estrai_serie_temporali``:
    NON ricalcola alcuna fisica, ricampiona (interpolazione lineare,
    ``numpy.interp``) le stesse serie gia' validate a istanti equispaziati
    nel tempo, aggiungendo il solo campo derivato ``angolo_volo_gradi``.

    Motivazione: i dati grezzi hanno solo ~57 punti totali (passo
    adattivo dell'integratore, molto irregolare tra le 3 fasi -- pochi
    punti in Fase A, molti nel gravity turn/tangente lineare).
    Sottocampionati per INDICE (come fa ``anima_traiettoria``), con cosi'
    pochi punti sorgente producono ripetizioni/scatti percettibili in
    un'animazione con un simbolo che deve ruotare visibilmente. Qui il
    ricampionamento e' nel TEMPO, non negli indici -- stesso principio di
    "e' rendering, non un nuovo risultato fisico" gia' dichiarato per la
    sottocampionatura per indice esistente.

    L'istante esatto di staging (seconda transizione di
    ``estrai_serie_temporali``) e' preservato come frame ESATTO (inserito
    esplicitamente nel vettore dei tempi campione, non solo approssimato
    dall'interpolazione), cosi' il cambio visivo allo staging cade
    esattamente nel frame giusto, non a un frame vicino per caso.

    Parametri
    ---------
    assemblaggio : dict
        Stesso input di ``estrai_serie_temporali``.
    n_frame : int, opzionale
        Numero di istanti equispaziati nel tempo (oltre all'istante di
        staging, sempre inserito esplicitamente). Default 220.

    Ritorna
    -------
    dict
        ``{"t", "x", "h", "v", "angolo_gradi", "idx_staging",
        "t_staging"}``: array numpy ricampionati a passo di tempo
        uniforme (lunghezza ``n_frame``, +1 se l'istante di staging non
        cadeva gia' su un campione). ``idx_staging`` e' l'indice, in
        questi array ricampionati, del frame che cade esattamente
        sull'istante di staging.
    """
    serie = estrai_serie_temporali(assemblaggio)
    t, x, h, v = serie["t"], serie["x"], serie["h"], serie["v"]

    angolo_gradi_raw = angolo_volo_gradi(assemblaggio, serie)
    t_staging = serie["t_transizioni"][1]

    # Tempi campione: n_frame equispaziati in [t[0], t[-1]], PIU' l'istante
    # esatto di staging inserito esplicitamente (evita che il cambio
    # visivo dipenda dalla risoluzione casuale del campionamento).
    t_campione = np.linspace(t[0], t[-1], n_frame)
    t_campione = np.sort(np.unique(np.concatenate([t_campione, [t_staging]])))

    x_i = np.interp(t_campione, t, x)
    h_i = np.interp(t_campione, t, h)
    v_i = np.interp(t_campione, t, v)
    angolo_i = np.interp(t_campione, t, angolo_gradi_raw)

    idx_staging = int(np.searchsorted(t_campione, t_staging))

    return {
        "t": t_campione,
        "x": x_i,
        "h": h_i,
        "v": v_i,
        "angolo_gradi": angolo_i,
        "idx_staging": idx_staging,
        "t_staging": t_staging,
    }


def grafico_quota_velocita_tempo(serie, percorso_output):
    """Sottografici h(t), |v|(t), m(t) con linee verticali alle transizioni di fase.

    Salva un PNG in ``percorso_output`` (crea la cartella se non esiste).

    Ritorna
    -------
    str
        ``percorso_output``, per comodita' di chaining/logging.
    """
    _assicura_cartella(percorso_output)

    fig, assi = plt.subplots(3, 1, sharex=True, figsize=(9, 10))

    assi[0].plot(serie["t"], serie["h"] / 1000.0, color="tab:blue")
    assi[0].set_ylabel("Quota h [km]")
    assi[0].set_title("Ascesa guidata: quota, velocita', massa nel tempo")

    assi[1].plot(serie["t"], serie["v"], color="tab:orange")
    assi[1].set_ylabel("Velocita' |v| [m/s]")

    assi[2].plot(serie["t"], serie["m"] / 1000.0, color="tab:green")
    assi[2].set_ylabel("Massa m [t]")
    assi[2].set_xlabel("Tempo t [s]")

    for ax in assi:
        for t_transizione, etichetta in zip(serie["t_transizioni"], serie["label_transizioni"]):
            ax.axvline(t_transizione, color="gray", linestyle="--", linewidth=1)
        ax.grid(alpha=0.3)

    # Etichette delle transizioni una sola volta, sul sottografico in alto.
    for t_transizione, etichetta in zip(serie["t_transizioni"], serie["label_transizioni"]):
        assi[0].annotate(
            etichetta,
            xy=(t_transizione, assi[0].get_ylim()[1]),
            xytext=(2, -10),
            textcoords="offset points",
            rotation=90,
            va="top",
            fontsize=8,
            color="gray",
        )

    fig.tight_layout()
    fig.savefig(percorso_output, dpi=150)
    plt.close(fig)
    return percorso_output


def grafico_traiettoria(serie, percorso_output):
    """Traiettoria nel piano (x, h). Salva un PNG in ``percorso_output``."""
    _assicura_cartella(percorso_output)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(serie["x"] / 1000.0, serie["h"] / 1000.0, color="tab:blue")

    for idx_transizione, etichetta in zip(serie["idx_transizioni"], serie["label_transizioni"]):
        ax.plot(
            serie["x"][idx_transizione] / 1000.0,
            serie["h"][idx_transizione] / 1000.0,
            "o",
            color="gray",
            label=etichetta,
        )

    ax.set_xlabel("Distanza a terra x [km]")
    ax.set_ylabel("Quota h [km]")
    ax.set_title("Traiettoria di ascesa nel piano verticale (x, h)")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(percorso_output, dpi=150)
    plt.close(fig)
    return percorso_output


def anima_traiettoria(serie, percorso_output, n_frame=150, fps=20):
    """Anima la traiettoria (x, h) nel tempo. Salva una GIF con ``PillowWriter``.

    Sottocampiona a ``n_frame`` frame equispaziati sugli INDICI disponibili
    (non sul tempo, che ha passo non uniforme per via del passo adattivo
    dell'integratore) per un'animazione fluida senza rallentamenti nei
    tratti con pochi punti sorgente (es. la Fase A, ~5 punti).

    Ritorna
    -------
    str
        ``percorso_output``.
    """
    _assicura_cartella(percorso_output)

    x_km = serie["x"] / 1000.0
    h_km = serie["h"] / 1000.0
    n_punti = len(serie["t"])

    indici_frame = np.linspace(0, n_punti - 1, min(n_frame, n_punti)).astype(int)

    fig, ax = plt.subplots(figsize=(8, 6))
    margine_x = max(x_km.max() * 0.05, 1.0)
    ax.set_xlim(x_km.min() - margine_x, x_km.max() + margine_x)
    ax.set_ylim(0, h_km.max() * 1.05)
    ax.set_xlabel("Distanza a terra x [km]")
    ax.set_ylabel("Quota h [km]")
    ax.set_title("Traiettoria di ascesa (animazione)")
    ax.grid(alpha=0.3)

    (linea,) = ax.plot([], [], "b-", linewidth=1.5)
    (punto,) = ax.plot([], [], "ro", markersize=6)

    def init():
        linea.set_data([], [])
        punto.set_data([], [])
        return linea, punto

    def aggiorna(frame_idx):
        i = indici_frame[frame_idx]
        linea.set_data(x_km[: i + 1], h_km[: i + 1])
        punto.set_data([x_km[i]], [h_km[i]])
        return linea, punto

    animazione = FuncAnimation(
        fig, aggiorna, init_func=init, frames=len(indici_frame), interval=1000.0 / fps, blit=True
    )
    animazione.save(percorso_output, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return percorso_output


# Forma del simbolo del razzo e della fiamma, in "spazio marker" (coordinate
# locali normalizzate, non unita' dati): punta lungo +x (orizzontale, verso
# destra) quando l'angolo di rotazione applicato e' 0. La rotazione a
# runtime (Affine2D().rotate_deg(angolo), angolo = angolo_volo_gradi) avviene
# in questo stesso spazio locale, PRIMA della trasformazione a coordinate
# schermo -- per questo il simbolo appare orientato correttamente anche con
# assi x/h a scale molto diverse (qui ~1593 km vs ~202 km): la rotazione non
# e' soggetta alla distorsione di aspect ratio dei dati.
_FORMA_RAZZO = Path(
    [(1.0, 0.0), (-0.6, 0.5), (-0.3, 0.0), (-0.6, -0.5), (1.0, 0.0)], closed=True
)
_FORMA_FIAMMA = Path(
    [(0.0, 0.0), (-1.0, 0.35), (-1.6, 0.0), (-1.0, -0.35), (0.0, 0.0)], closed=True
)


def anima_traiettoria_avanzata(serie_interpolata, percorso_output, fps=30):
    """Animazione stilizzata della traiettoria: razzo orientato, scia, stage separation.

    Si AFFIANCA a ``anima_traiettoria`` (non la sostituisce, non la
    modifica): stessa filosofia "nessuna nuova fisica, solo rendering di
    numeri gia' calcolati e validati", con un simbolo del razzo che ruota
    visibilmente seguendo l'angolo di volo reale (``angolo_volo_gradi``),
    una scia con dissolvenza che mostra la traiettoria gia' percorsa, e un
    cambio di colore/dimensione preciso all'istante esatto di stage
    separation (``serie_interpolata["idx_staging"]``, un frame REALE del
    ricampionamento, non approssimato).

    Il simbolo del razzo e della fiamma sono disegnati a dimensione FISSA
    IN PUNTI SCHERMO (``markersize``, non unita' dati): un razzo alla sua
    scala reale (~70 m) sarebbe un punto invisibile su un asse di
    centinaia di km. E' una scelta di resa visiva, non un dato fisico --
    dichiarata esplicitamente qui per lo stesso principio di onesta' sulle
    approssimazioni gia' seguito nel resto del progetto (li' per la fisica,
    qui per il rendering). La rotazione del simbolo avviene in spazio
    marker (vedi ``_FORMA_RAZZO``/``_FORMA_FIAMMA``), quindi resta corretta
    a schermo indipendentemente dall'aspect ratio degli assi.

    Parametri
    ---------
    serie_interpolata : dict
        Risultato di ``estrai_serie_temporali_interpolata``.
    percorso_output : str
        Percorso del file GIF di output.
    fps : int, opzionale
        Fotogrammi al secondo della GIF di output. Default 30.

    Ritorna
    -------
    str
        ``percorso_output``.
    """
    _assicura_cartella(percorso_output)

    t = serie_interpolata["t"]
    x_km = serie_interpolata["x"] / 1000.0
    h_km = serie_interpolata["h"] / 1000.0
    v = serie_interpolata["v"]
    angolo = serie_interpolata["angolo_gradi"]
    idx_staging = serie_interpolata["idx_staging"]
    n = len(t)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    margine_x = max((x_km.max() - x_km.min()) * 0.05, 1.0)
    margine_h = max(h_km.max() * 0.08, 1.0)
    lim_x = (x_km.min() - margine_x, x_km.max() + margine_x)
    lim_h = (-margine_h * 0.3, h_km.max() + margine_h)
    ax.set_xlim(*lim_x)
    ax.set_ylim(*lim_h)
    ax.set_xlabel("Distanza a terra x [km]")
    ax.set_ylabel("Quota h [km]")
    ax.set_title("Ascesa guidata — animazione")

    # --- Sfondo: gradiente atmosfera -> spazio, puramente estetico. ---
    # Tappe scelte in ordine di luminosita' STRETTAMENTE decrescente (ogni
    # canale RGB decresce ad ogni tappa) per un imbrunimento monotono dal
    # blu atmosferico vicino al suolo al nero quasi assoluto in quota --
    # un gradiente non monotono (es. una tappa piu' chiara in mezzo)
    # produrrebbe una fascia illuminata artificiale a meta' altezza.
    gradiente = np.linspace(0, 1, 256).reshape(-1, 1)
    cmap_cielo = matplotlib.colors.LinearSegmentedColormap.from_list(
        "cielo_notturno", ["#1b3a63", "#0f2444", "#081326", "#02040a"]
    )
    ax.imshow(
        gradiente,
        extent=[lim_x[0], lim_x[1], lim_h[0], lim_h[1]],
        aspect="auto",
        cmap=cmap_cielo,
        zorder=0,
        origin="lower",
    )

    # --- Terreno. ---
    ax.axhspan(lim_h[0], 0, color="#2a1f14", zorder=1)
    ax.axhline(0, color="#5c4326", linewidth=1.5, zorder=1)
    ax.grid(alpha=0.15, color="white", zorder=1)

    colore_stadio1 = "#4fc3f7"
    colore_stadio2 = "#ff8a3d"

    scia = LineCollection([], zorder=2)
    ax.add_collection(scia)

    (fiamma,) = ax.plot(
        [],
        [],
        marker=MarkerStyle(_FORMA_FIAMMA),
        markersize=20,
        markerfacecolor="#ffb347",
        markeredgewidth=0,
        linestyle="",
        alpha=0.9,
        zorder=3,
    )
    (razzo,) = ax.plot(
        [],
        [],
        marker=MarkerStyle(_FORMA_RAZZO),
        markersize=15,
        markerfacecolor="white",
        markeredgecolor="#333333",
        markeredgewidth=0.8,
        linestyle="",
        zorder=4,
    )

    testo_flash = ax.text(
        0.5,
        0.93,
        "",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=13,
        color="#ffdd55",
        fontweight="bold",
        zorder=5,
    )
    hud = ax.text(
        0.02,
        0.96,
        "",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="white",
        family="monospace",
        zorder=5,
        bbox=dict(boxstyle="round", facecolor="black", alpha=0.45, edgecolor="none"),
    )

    # Ampiezza della finestra di "flash" attorno allo stage separation
    # (in numero di frame, non in secondi: dipende da fps solo tramite il
    # numero di frame gia' fissato da estrai_serie_temporali_interpolata).
    finestra_flash = max(int(0.4 * fps), 3)

    def init():
        scia.set_segments([])
        razzo.set_data([], [])
        fiamma.set_data([], [])
        testo_flash.set_text("")
        hud.set_text("")
        return scia, razzo, fiamma, testo_flash, hud

    def aggiorna(i):
        # Scia con dissolvenza: alpha crescente verso il presente (effetto
        # "cometa"), colorata per stadio (cambio di colore ESATTO
        # all'istante di staging, non un'approssimazione).
        if i >= 1:
            punti = np.column_stack([x_km[: i + 1], h_km[: i + 1]])
            segmenti = np.stack([punti[:-1], punti[1:]], axis=1)
            n_seg = len(segmenti)
            alphas = np.linspace(0.15, 0.95, n_seg) if n_seg > 1 else [0.95]
            colori_rgba = []
            for k in range(n_seg):
                colore = colore_stadio1 if k < idx_staging else colore_stadio2
                r, g, b = matplotlib.colors.to_rgb(colore)
                colori_rgba.append((r, g, b, alphas[k]))
            scia.set_segments(segmenti)
            scia.set_color(colori_rgba)
            scia.set_linewidth(3.0)
        else:
            scia.set_segments([])

        # Razzo/fiamma orientati lungo l'angolo di volo reale: rotazione in
        # spazio marker (screen-space), non in unita' dati (vedi docstring
        # e commento su _FORMA_RAZZO/_FORMA_FIAMMA).
        theta = angolo[i]
        razzo.set_marker(MarkerStyle(_FORMA_RAZZO, transform=Affine2D().rotate_deg(theta)))
        razzo.set_data([x_km[i]], [h_km[i]])
        fiamma.set_marker(MarkerStyle(_FORMA_FIAMMA, transform=Affine2D().rotate_deg(theta)))
        fiamma.set_data([x_km[i]], [h_km[i]])

        # Flash allo stage separation: simbolo temporaneamente piu' grande
        # e di colore diverso per una breve finestra di frame CENTRATA
        # sul frame esatto di staging (idx_staging, non approssimato).
        if abs(i - idx_staging) <= finestra_flash:
            razzo.set_markersize(24)
            razzo.set_markerfacecolor("#ffdd55")
            testo_flash.set_text("STAGE SEPARATION")
        else:
            razzo.set_markersize(15)
            razzo.set_markerfacecolor("white")
            testo_flash.set_text("")

        stadio_corrente = 1 if i < idx_staging else 2
        hud.set_text(
            f"t = {t[i]:6.1f} s\nh = {h_km[i]:7.2f} km\nv = {v[i]:7.1f} m/s\nstadio {stadio_corrente}"
        )

        return scia, razzo, fiamma, testo_flash, hud

    animazione = FuncAnimation(fig, aggiorna, init_func=init, frames=n, interval=1000.0 / fps, blit=False)
    animazione.save(percorso_output, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return percorso_output


def riepilogo_testuale(risultato_validazione):
    """Riepilogo numerico testuale, riusa SENZA RICALCOLO i numeri di ``esegui_validazione``.

    Nessuna formula di delta-v/perdite reimplementata qui: legge
    direttamente i campi gia' calcolati da
    ``lanciatore.validazione.calcola_delta_v_ideale`` /
    ``lanciatore.validazione.scompone_perdite`` (via
    ``esegui_validazione``), esattamente come richiesto dal piano
    (STATUS.md, Ciclo 9).

    Parametri
    ---------
    risultato_validazione : dict
        Risultato di ``lanciatore.validazione.esegui_validazione()``.

    Ritorna
    -------
    str
        Riepilogo multi-riga, pronto per la stampa.
    """
    dv = risultato_validazione["delta_v_ideale"]
    perdite = risultato_validazione["perdite"]
    assemblaggio = risultato_validazione["assemblaggio"]

    righe = [
        "Riepilogo numerico (Step 7, valori riusati senza ricalcolo):",
        f"  delta-v ideale totale (Tsiolkovsky): {dv['dv_totale']:.1f} m/s"
        f" (stadio 1: {dv['dv1']:.1f} m/s, stadio 2: {dv['dv2']:.1f} m/s)",
        f"  velocita' finale raggiunta: {assemblaggio['v_finale']:.1f} m/s"
        f" (target orizzontale: {assemblaggio['vx_target']:.1f} m/s,"
        f" residuo: {assemblaggio['residuo_modulo']:.1f} m/s)",
        "  scomposizione delle perdite:"
        f" gravita'={perdite['perdita_gravita']:.1f} m/s,"
        f" drag={perdite['perdita_drag']:.1f} m/s,"
        f" manovra={perdite['perdita_manovra']:.1f} m/s",
        f"  residuo di chiusura dell'identita': {perdite['residuo_chiusura']:.3e} m/s"
        f" ({perdite['residuo_chiusura_relativo']:.2e} relativo)",
    ]
    return "\n".join(righe)


if __name__ == "__main__":
    # Genera effettivamente i 3 file di output sui dati di
    # esegui_validazione() (STATUS.md Ciclo 9, task 4) e stampa il
    # riepilogo numerico testuale.
    risultato = val.esegui_validazione()
    serie = estrai_serie_temporali(risultato["assemblaggio"])
    serie_interpolata = estrai_serie_temporali_interpolata(risultato["assemblaggio"])

    percorso_grafico_tempo = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "quota_velocita_tempo.png")
    percorso_grafico_traiettoria = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "traiettoria.png")
    percorso_animazione = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "traiettoria.gif")
    percorso_animazione_avanzata = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "traiettoria_avanzata.gif")

    grafico_quota_velocita_tempo(serie, percorso_grafico_tempo)
    grafico_traiettoria(serie, percorso_grafico_traiettoria)
    anima_traiettoria(serie, percorso_animazione)
    anima_traiettoria_avanzata(serie_interpolata, percorso_animazione_avanzata)

    print(riepilogo_testuale(risultato))
    print()
    print("File generati:")
    print(f"  {percorso_grafico_tempo}")
    print(f"  {percorso_grafico_traiettoria}")
    print(f"  {percorso_animazione}")
    print(f"  {percorso_animazione_avanzata}")
