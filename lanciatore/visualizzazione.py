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
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

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

    percorso_grafico_tempo = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "quota_velocita_tempo.png")
    percorso_grafico_traiettoria = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "traiettoria.png")
    percorso_animazione = os.path.join(_CARTELLA_OUTPUT_PROGETTO, "traiettoria.gif")

    grafico_quota_velocita_tempo(serie, percorso_grafico_tempo)
    grafico_traiettoria(serie, percorso_grafico_traiettoria)
    anima_traiettoria(serie, percorso_animazione)

    print(riepilogo_testuale(risultato))
    print()
    print("File generati:")
    print(f"  {percorso_grafico_tempo}")
    print(f"  {percorso_grafico_traiettoria}")
    print(f"  {percorso_animazione}")
