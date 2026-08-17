"""Staging durante la fase di guida a tangente lineare (Step 8).

Scenario dimostrativo (vedi STATUS.md, Ciclo 8 + addendum reviewer): un
Segmento 1 pianifica ``A1, B1`` per un target di velocita' finale
assumendo (erroneamente) di bruciare fino a un istante nominale
``TF1_NOMINALE``, ma il propellente realmente disponibile a bordo
esaurisce PRIMA di quell'istante — vincolo fisico reale del veicolo, non
un taglio arbitrario. Lo staging e' rilevato come un vero evento
terminale dell'integrazione ODE (stesso pattern di
``evento_fine_propellente_2d`` in ``lanciatore/staging.py``), non un
cutoff scriptato a un tempo scelto a mano. Il Segmento 2 deve ricalcolare
``A2, B2`` dallo stato EFFETTIVO di handoff per raggiungere comunque lo
stesso target: questo modulo dimostra sia perche' il ricalcolo e'
necessario (i vecchi coefficienti, continuati a tempo assoluto invariato,
mancano il target) sia che il ricalcolo funziona.

Riuso totale di ``lanciatore/guida_esoatmosferica.py``
(``derivate_stato_tangente_lineare``, ``integra_tangente_lineare``,
``risolvi_coefficienti_tangente_lineare``), nessuna modifica a quel
modulo. ``lanciatore/staging.py`` NON viene importato: e' usato solo come
riferimento di stile per il pattern degli eventi ODE (evento terminale
genuino via ``functools.partial`` sulla soglia di massa, stesso principio
di continuita' x/h/vx/vh "riporto diretto" gia' documentato li' per lo
staging in fase atmosferica).

Parametri dello scenario (Ciclo 8, verificati numericamente prima del
piano — vedi STATUS.md)
---------------------------------------------------------------------
- Segmento 1: ``m0=1000`` kg, ``spinta=40000`` N, ``isp=300`` s,
  ``h0=5000`` m (margine quota positivo, stesso principio di Fase A del
  Ciclo 7). Target finale ``vx=1800`` m/s, ``vh=0`` m/s (entro capacita'
  ideale Tsiolkovsky del solo Segmento 1). ``tf1`` nominale usato per la
  PIANIFICAZIONE di ``A1, B1`` = 40 s (ipotesi di bruciamento completo,
  poi rivelatasi errata).
- Propellente REALMENTE disponibile nel Segmento 1: solo 15 s di
  bruciamento (``mdot1 = spinta/(isp*G0)``), quindi
  ``m_prop1_reale = mdot1*15.0``, ``m_soglia1 = m0 - m_prop1_reale``.
  L'integrazione del Segmento 1 si ferma qui per un evento terminale
  genuino, non a ``t=15`` per costruzione scriptata: 15 s e' semplicemente
  il tempo che la fisica (portata massica costante) impiega a raggiungere
  quella soglia di massa, verificato a posteriori tramite ``t_events``,
  non imposto come limite di integrazione.
- Discontinuita' di massa allo staging: massa EFFETTIVA restituita da
  ``solve_ivp`` a fine Segmento 1 meno ``m_strut1=100.0`` kg (provvisorio)
  = massa di ignizione del Segmento 2. ``x, h, vx, vh`` sono riportati
  DIRETTAMENTE (nessuna conversione: qui, a differenza dello Step 5/7, lo
  stato e' gia' cartesiano su entrambi i lati dello staging).
- Segmento 2: ``spinta=40000`` N, ``isp=320`` s, tempo residuo
  ``tf2=25`` s. I coefficienti ``A2, B2`` sono ricalcolati dallo stato di
  handoff con ``risolvi_coefficienti_tangente_lineare``, usando come
  guess ``A0=A1, B0=B1`` — regola di seeding DETERMINISTICA (non un
  guess arbitrario): i coefficienti convergenti del segmento precedente
  sono il punto di partenza del root-finder per il segmento successivo,
  stesso principio gia' scritto nella docstring di
  ``risolvi_coefficienti_tangente_lineare`` ("il PEG non riparte mai da
  un default generico, riusa la soluzione del ciclo precedente").

Prova di invalidita' dei vecchi coefficienti (criterio 2 del piano)
---------------------------------------------------------------------
"Cosa sarebbe successo SENZA ricalcolo" significa continuare il
MEDESIMO profilo di guida ``theta(t) = arctan(A1 + B1*t)`` a tempo
ASSOLUTO invariato (l'orologio della legge di guida non riparte da zero
al Segmento 2) — equivalente a integrare il Segmento 2 con
``A0_continuato = A1 + B1*t_stage`` (traslazione dell'offset temporale
dello staging dentro la costante additiva) e ``B0_continuato = B1``
(stessa pendenza). Questo NON e' lo stesso calcolo di "integrare con A1,
B1 invariati e t che riparte da zero" (quello sarebbe un artefatto di
reset del tempo, discusso e scartato nell'addendum del reviewer, Ciclo 8):
qui si usa deliberatamente la continuazione a tempo assoluto perche' e'
quella fisicamente corretta (lo stesso profilo di guida pianificato in
origine, semplicemente eseguito con lo stato di handoff reale invece di
quello ipotizzato in fase di pianificazione).

Fenomeno delle radici multiple in questo contesto (vedi anche
``risolvi_coefficienti_tangente_lineare``, Step 6): la regola di seeding
``A0=A1, B0=B1`` converge in modo affidabile alla radice fisicamente
continua. Un guess diverso e non imparentato con il segmento precedente
(es. ``A0=0.5, B0=-0.01``) puo' convergere a una radice DIVERSA ma
ugualmente valida (stesso target raggiunto, coefficienti diversi) — non
e' un bug del risolutore, e' lo stesso fenomeno gia' documentato allo
Step 6, qui verificato di nuovo in un contesto a due segmenti (vedi
``tests/test_staging_esoatmosferico.py``).
"""

import functools

import numpy as np
from scipy.integrate import solve_ivp

from lanciatore.costanti import G0
from lanciatore.guida_esoatmosferica import (
    derivate_stato_tangente_lineare,
    integra_tangente_lineare,
    risolvi_coefficienti_tangente_lineare,
)

# ---------------------------------------------------------------------------
# Parametri dello scenario dimostrativo (Ciclo 8), verificati numericamente
# nel piano prima dell'implementazione. Non sono costanti fisiche di
# progetto (restano qui, non in lanciatore/costanti.py), stesso principio
# gia' seguito per i parametri veicolo "provvisori" degli step precedenti
# (es. Step 2, Step 5).
# ---------------------------------------------------------------------------
M0_SEG1 = 1000.0  # kg, massa iniziale Segmento 1
SPINTA_SEG1 = 40000.0  # N
ISP_SEG1 = 300.0  # s
H0_SEG1 = 5000.0  # m, margine di quota positivo (stesso principio di Fase A, Ciclo 7)

VX_TARGET = 1800.0  # m/s
VH_TARGET = 0.0  # m/s

TF1_NOMINALE = 40.0  # s, orizzonte assunto (erroneamente) in fase di pianificazione
A0_GUESS_SEG1 = 0.2  # guess iniziale per A1 (verificato: converge)
B0_GUESS_SEG1 = 0.0  # guess iniziale per B1 (verificato: converge)

DURATA_BRUCIAMENTO_REALE_SEG1 = 15.0  # s, propellente realmente a bordo (vincolo fisico)
M_STRUT1 = 100.0  # kg, massa strutturale espulsa allo staging (provvisorio)

SPINTA_SEG2 = 40000.0  # N
ISP_SEG2 = 320.0  # s
TF2 = 25.0  # s, tempo residuo nominale per il Segmento 2 (40 - 15)


def _mdot(spinta, isp):
    """Portata massica costante, ``spinta/(isp*G0)`` (convenzione internazionale, G0 costante)."""
    return spinta / (isp * G0)


def _evento_fine_propellente_reale(t, y, A, B, spinta, mdot, m_soglia):
    """Evento terminale genuino: si annulla quando la massa scende a ``m_soglia``.

    Stessa firma di ``derivate_stato_tangente_lineare`` (t, y, A, B,
    spinta, mdot) piu' ``m_soglia`` legato con ``functools.partial`` prima
    di passare l'evento a ``solve_ivp`` — stesso pattern di
    ``evento_fine_propellente_2d`` in ``lanciatore/staging.py``: ``args``
    e' condiviso tra la funzione derivate e le funzioni evento passate a
    ``solve_ivp``, quindi l'evento deve accettare gli stessi argomenti
    posizionali, con l'argomento extra legato a monte.
    """
    return y[4] - m_soglia


def esegui_staging_esoatmosferico(
    m0_seg1=M0_SEG1,
    spinta_seg1=SPINTA_SEG1,
    isp_seg1=ISP_SEG1,
    h0_seg1=H0_SEG1,
    vx_target=VX_TARGET,
    vh_target=VH_TARGET,
    tf1_nominale=TF1_NOMINALE,
    a0_guess_seg1=A0_GUESS_SEG1,
    b0_guess_seg1=B0_GUESS_SEG1,
    durata_bruciamento_reale_seg1=DURATA_BRUCIAMENTO_REALE_SEG1,
    m_strut1=M_STRUT1,
    spinta_seg2=SPINTA_SEG2,
    isp_seg2=ISP_SEG2,
    tf2=TF2,
):
    """Esegue lo scenario di staging in fase esoatmosferica (Step 8).

    Segue esattamente il design del piano consolidato (STATUS.md, Ciclo
    8 + addendum reviewer): pianifica A1,B1 per il target finale
    assumendo un bruciamento completo fino a ``tf1_nominale``; integra il
    Segmento 1 fino all'evento terminale GENUINO di fine propellente
    reale (soglia di massa raggiunta a ``durata_bruciamento_reale_seg1``
    secondi di bruciamento a portata costante, non un cutoff scriptato);
    applica la discontinuita' di massa allo staging (massa effettiva meno
    ``m_strut1``); ricalcola A2,B2 dallo stato EFFETTIVO di handoff con
    guess deterministico ``A0=A1, B0=B1``; calcola anche il risultato
    "senza ricalcolo" (continuazione del piano originale a tempo
    assoluto) per dimostrare l'errore che il ricalcolo evita.

    Nessun parametro ha un default nascosto rispetto al piano: tutti i
    default corrispondono esattamente ai valori dello scenario
    verificato in STATUS.md, ma sono argomenti espliciti (non costanti
    module-level usate direttamente nel corpo della funzione) cosi' i
    test possono variarli per esplorare, ad esempio, la robustezza del
    fenomeno delle radici multiple senza duplicare la logica.

    Ritorna
    -------
    dict
        Nessun dato nascosto/appiattito:

        - ``risultato_seg1``: oggetto grezzo ``solve_ivp`` (OdeResult)
          del Segmento 1, integrato fino all'evento di fine propellente
          reale. Accesso diretto a ``.t``, ``.y``, ``.status``,
          ``.t_events`` per verificare la genuinita' dell'evento.
        - ``risultato_seg2``: oggetto grezzo ``solve_ivp`` del Segmento 2
          RICALCOLATO (A2, B2), che raggiunge il target.
        - ``risultato_senza_ricalcolo``: oggetto grezzo ``solve_ivp`` del
          Segmento 2 integrato con i coefficienti del piano ORIGINALE
          continuato a tempo assoluto (A0_continuato, B0_continuato),
          per il confronto esplicito che dimostra l'errore.
        - ``A1, B1``: coefficienti nominali del Segmento 1.
        - ``A2, B2``: coefficienti ricalcolati del Segmento 2.
        - ``A0_continuato, B0_continuato``: coefficienti "senza
          ricalcolo" usati per ``risultato_senza_ricalcolo``.
        - ``t_stage``: istante (tempo assoluto, s) in cui l'evento di
          fine propellente reale ha fermato il Segmento 1 (letto da
          ``risultato_seg1.t[-1]`` / ``t_events``, non hardcodato).
        - ``m_eff``: massa effettiva (kg) restituita da ``solve_ivp`` a
          fine Segmento 1 (prima della sottrazione della massa
          strutturale).
        - ``m_strut1``: massa strutturale espulsa allo staging, kg.
        - ``m_ignizione2``: massa di ignizione del Segmento 2
          (``m_eff - m_strut1``), kg.
        - ``x_eff, h_eff, vx_eff, vh_eff``: stato di posizione/velocita'
          allo staging (identico, per riporto diretto, tra fine
          Segmento 1 e inizio Segmento 2).
        - ``vx_target, vh_target``: target di velocita' (ripassati per
          comodita' del chiamante/test, stessi valori degli argomenti).
        - ``mdot1, mdot2``: portate massiche costanti dei due segmenti,
          kg/s.
        - ``m_ignizione2, spinta_seg2, mdot2, tf2, x_eff, h_eff, vx_eff,
          vh_eff``: ripassati esplicitamente anche come chiavi singole
          (oltre che deducibili da quanto sopra) per permettere ai test
          di richiamare ``risolvi_coefficienti_tangente_lineare`` con un
          guess diverso senza dover ricostruire lo stato di handoff a
          mano.
    """
    mdot1 = _mdot(spinta_seg1, isp_seg1)

    # --- Segmento 1: pianificazione A1,B1 assumendo (erroneamente) il
    # bruciamento completo fino a tf1_nominale. ---
    A1, B1 = risolvi_coefficienti_tangente_lineare(
        m0_seg1,
        spinta_seg1,
        mdot1,
        tf1_nominale,
        vx_target,
        vh_target,
        a0_guess_seg1,
        b0_guess_seg1,
        h0=h0_seg1,
    )

    # --- Integrazione del Segmento 1 fino all'evento GENUINO di fine
    # propellente reale (soglia di massa raggiunta fisicamente prima di
    # tf1_nominale, non un t_span troncato a mano). solve_ivp DIRETTO
    # (non integra_tangente_lineare, che non espone `events`). ---
    m_prop1_reale = mdot1 * durata_bruciamento_reale_seg1
    m_soglia1 = m0_seg1 - m_prop1_reale

    evento_seg1 = functools.partial(
        _evento_fine_propellente_reale, m_soglia=m_soglia1
    )
    evento_seg1.terminal = True
    evento_seg1.direction = -1

    y0_seg1 = np.array([0.0, h0_seg1, 0.0, 0.0, m0_seg1])
    risultato_seg1 = solve_ivp(
        derivate_stato_tangente_lineare,
        t_span=(0.0, tf1_nominale),
        y0=y0_seg1,
        method="RK45",
        args=(A1, B1, spinta_seg1, mdot1),
        events=evento_seg1,
        rtol=1e-8,
        # Stesse scale/tolleranze di integra_tangente_lineare (Step 4/6).
        atol=[1e-3, 1e-3, 1e-6, 1e-6, 1e-3],
        dense_output=False,
    )

    x_eff, h_eff, vx_eff, vh_eff, m_eff = risultato_seg1.y[:, -1]
    t_stage = float(risultato_seg1.t[-1])

    # --- Discontinuita' di massa allo staging: riporto diretto di
    # x,h,vx,vh; massa di ignizione del Segmento 2 = massa effettiva
    # meno la massa strutturale espulsa. ---
    m_ignizione2 = m_eff - m_strut1

    # --- Segmento 2: ricalcolo di A2,B2 dallo stato EFFETTIVO di
    # handoff, guess deterministico A0=A1, B0=B1 (regola di seeding,
    # non un guess arbitrario). ---
    mdot2 = _mdot(spinta_seg2, isp_seg2)
    A2, B2 = risolvi_coefficienti_tangente_lineare(
        m_ignizione2,
        spinta_seg2,
        mdot2,
        tf2,
        vx_target,
        vh_target,
        A1,
        B1,
        x0=x_eff,
        h0=h_eff,
        vx0=vx_eff,
        vh0=vh_eff,
    )

    risultato_seg2 = integra_tangente_lineare(
        m_ignizione2,
        spinta_seg2,
        mdot2,
        A2,
        B2,
        tf2,
        x0=x_eff,
        h0=h_eff,
        vx0=vx_eff,
        vh0=vh_eff,
    )

    # --- "Cosa sarebbe successo SENZA ricalcolo": continuazione del
    # piano ORIGINALE a tempo ASSOLUTO (l'orologio della legge di guida
    # non riparte da zero al Segmento 2). ---
    A0_continuato = A1 + B1 * t_stage
    B0_continuato = B1
    risultato_senza_ricalcolo = integra_tangente_lineare(
        m_ignizione2,
        spinta_seg2,
        mdot2,
        A0_continuato,
        B0_continuato,
        tf2,
        x0=x_eff,
        h0=h_eff,
        vx0=vx_eff,
        vh0=vh_eff,
    )

    return {
        "risultato_seg1": risultato_seg1,
        "risultato_seg2": risultato_seg2,
        "risultato_senza_ricalcolo": risultato_senza_ricalcolo,
        "A1": A1,
        "B1": B1,
        "A2": A2,
        "B2": B2,
        "A0_continuato": A0_continuato,
        "B0_continuato": B0_continuato,
        "t_stage": t_stage,
        "m_eff": m_eff,
        "m_strut1": m_strut1,
        "m_ignizione2": m_ignizione2,
        "x_eff": x_eff,
        "h_eff": h_eff,
        "vx_eff": vx_eff,
        "vh_eff": vh_eff,
        "vx_target": vx_target,
        "vh_target": vh_target,
        "mdot1": mdot1,
        "mdot2": mdot2,
        "spinta_seg2": spinta_seg2,
        "tf2": tf2,
    }
