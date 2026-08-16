"""Multistadio, fase atmosferica: eventi di staging, cambio massa discontinuo (Step 5).

Scope (decisione utente 2026-08-16, vedi STATUS.md Ciclo 5): lo staging
qui riguarda SOLO la fase atmosferica (ascesa verticale + gravity turn,
Step 2/3), non la guida a tangente lineare (Step 4) — quella estensione
e' rimandata allo Step 7. Il meccanismo e' pero' implementato in forma
GENERICA e riusabile, cosi' la stessa struttura potra' servire anche
allo Step 7.

Rappresentazione di uno stadio: dict con ``m_prop`` (kg, propellente
bruciato dallo stadio), ``m_strut`` (kg, massa strutturale espulsa al
burnout), ``spinta`` (N), ``isp`` (s). Nessun campo "payload" esplicito:
resta implicito in ``m0`` totale meno tutto cio' che viene bruciato/
espulso dagli stadi elencati.

Riuso massimo, nessuna reimplementazione della fisica (vedi STATUS.md,
Ciclo 5, punto 1 "Design tecnico"):
    - Stadio 1: riusa DIRETTAMENTE ``guida.integra_gravity_turn`` (Fase
      A verticale + kick + Fase B gravity turn, Step 3), passando come
      ``m_vuoto`` la soglia di burnout DI QUESTO STADIO
      (``m0 - stadi[0]['m_prop']``): la funzione di Step 3 non ha
      bisogno di sapere che ci sono altri stadi attaccati sopra.
    - Stadi 2..N: riusano DIRETTAMENTE
      ``guida.derivate_stato_gravity_turn`` e i suoi eventi difensivi
      (``evento_impatto_suolo_2d``, ``evento_traiettoria_invalida``,
      ``evento_velocita_minima``), con un nuovo
      ``evento_fine_propellente_2d`` parametrizzato via
      ``functools.partial`` sulla soglia di burnout di QUEL segmento
      (stesso pattern gia' stabilito in ``guida.integra_gravity_turn``
      per Fase B).
    - Nessuna modifica a ``lanciatore/guida.py`` (diff vuoto).

Discontinuita' allo staging (vedi STATUS.md, Ciclo 5, addendum reviewer
punto 3): al burnout dello stadio i, ``x, h, v, gamma`` restano
ESATTAMENTE continui (presi dallo STATO EFFETTIVO restituito da
``solve_ivp``, ``.y[:, -1]``, non da un ricalcolo indipendente); la
massa di partenza dello stadio i+1 e' la massa EFFETTIVA al burnout
MENO ``m_strut`` dello stadio appena concluso — non il valore nominale
``m_vuoto_i`` calcolato a monte, per lo stesso principio gia' stabilito
per il kick di Step 3 (evitare un disallineamento con la tolleranza del
root-finder di ``solve_ivp``).

Nota su un evento di burnout intercettato da un evento DIFENSIVO: il
piano (STATUS.md, Ciclo 5, punto 3, criterio 4) affida esplicitamente
alla suite di test il compito di segnalare questo caso ("se succede e'
un segnale da investigare, non da ignorare"), non a un controllo che
interrompe l'esecuzione dentro il modulo stesso — questa funzione
percio' NON solleva errori se un segmento termina per un evento
difensivo (impatto suolo, traiettoria invalida, velocita' minima)
invece che per fine-propellente: prosegue con lo stato effettivo
comunque restituito da ``solve_ivp`` (che resta un numero finito, dato
che gli eventi sono ``terminal=True``), lasciando al chiamante/test la
diagnosi tramite ``t_events``. Aggiungere qui un controllo bloccante
sarebbe stata un'estensione non richiesta dal piano (vedi report dello
Step 5 in STATUS.md per la scoperta empirica che ha motivato questa
scelta).
"""

import functools

import numpy as np
from scipy.integrate import solve_ivp

from lanciatore import guida
from lanciatore.costanti import G0

# Ordine degli eventi nella tupla ``events`` passata a solve_ivp sia in
# guida.integra_gravity_turn (Fase B) sia nella chiamata diretta qui
# sotto per gli stadi 2..N: identico in entrambi i casi, usato da
# _evento_terminale per etichettare quale evento ha fermato ogni
# segmento.
_NOMI_EVENTI = (
    "fine_propellente",
    "impatto_suolo",
    "traiettoria_invalida",
    "velocita_minima",
)


def _evento_terminale(risultato):
    """Ritorna ``(nome_evento, istante)`` dell'evento che ha fermato il segmento.

    Usato per popolare il riepilogo in modo robusto anche nel caso (non
    nominale, ma da riportare esplicitamente, non da far fallire con un
    IndexError) in cui il segmento si sia fermato per un evento
    DIFENSIVO invece che per fine-propellente — vedi nota di modulo:
    la diagnosi di questo caso spetta al chiamante/test (criterio 4 del
    piano), non a un controllo bloccante qui.
    """
    for nome, t_ev in zip(_NOMI_EVENTI, risultato.t_events):
        if len(t_ev) > 0:
            return nome, t_ev[0]
    return None, risultato.t[-1]


def integra_multistadio_gravity_turn(
    m0,
    stadi,
    cd,
    area,
    v_kick,
    kick_angle_deg,
    h0=0.0,
    t_max=1000.0,
    v_min=1.0,
):
    """Integra l'ascesa multistadio in fase atmosferica (verticale + gravity turn).

    Parametri
    ---------
    m0 : float
        Massa totale iniziale allo stacco, kg (tutti gli stadi + payload
        implicito).
    stadi : list of dict
        Uno per stadio, nell'ordine di accensione, ciascuno con le
        chiavi ``m_prop`` (kg), ``m_strut`` (kg), ``spinta`` (N),
        ``isp`` (s).
    cd, area : float
        Coefficiente di resistenza aerodinamica e area di riferimento,
        GLOBALI per l'intera simulazione (decisione di design esplicita,
        addendum reviewer Ciclo 5 punto 2: coerente con "Cd costante"
        gia' dichiarato a livello di progetto in CLAUDE.md — differenze
        aerodinamiche tra stadi, es. sgancio ogiva, sono un
        raffinamento esplicitamente fuori scope qui).
    v_kick, kick_angle_deg : float
        Parametri della legge di guida (kick angle) usati per lo stadio
        1, stesso significato di ``guida.integra_gravity_turn``.
    h0 : float, opzionale
        Quota iniziale, m (default 0.0).
    t_max : float, opzionale
        Tempo massimo di integrazione PER CIASCUN SEGMENTO, s (default
        1000.0): limite di sicurezza, ogni segmento si ferma prima per
        uno dei suoi eventi terminali in condizioni normali.
    v_min : float, opzionale
        Soglia di velocita' minima per l'evento difensivo
        ``evento_velocita_minima`` in ogni segmento di gravity turn,
        m/s (default 1.0).

    Ritorna
    -------
    dict
        ``{"risultati": risultati, "riepilogo": riepilogo}``.

        ``risultati`` e' una lista con UN ELEMENTO PER STADIO (nessun
        risultato nascosto): l'elemento dello stadio 1 e' il dict grezzo
        restituito da ``guida.integra_gravity_turn`` (contiene sia
        ``"fase_verticale"`` sia ``"gravity_turn"``, cosi' com'e', senza
        appiattirlo); gli elementi degli stadi successivi sono
        direttamente gli oggetti risultato grezzi di ``solve_ivp`` per
        quel segmento di gravity turn (nessuna Fase A per gli stadi
        2..N: lo staging avviene sempre a velocita' ben superiore a
        ``v_kick``, quindi non serve e non ha senso fisico ripetere
        l'ascesa verticale pura).

        ``riepilogo`` e' un dict con una entry ``"stadio_i"`` (i da 1 a
        N) per ciascuno stadio, contenente ``massa_burnout_effettiva``,
        ``massa_burnout_nominale`` (soglia ``m_vuoto`` nominale di quel
        segmento), ``differenza_massa`` (effettiva meno nominale) e
        ``t_burnout`` (istante, nel tempo LOCALE di quel segmento,
        dell'evento di fine propellente).

    Solleva
    -------
    Nessuna eccezione specifica di questa funzione. Se un segmento
    termina per un evento difensivo (impatto suolo, traiettoria
    invalida, velocita' minima) invece che per fine propellente,
    l'integrazione prosegue comunque con lo stato effettivo restituito
    da ``solve_ivp`` (finito per costruzione): la diagnosi di questo
    caso spetta al chiamante/test tramite ``riepilogo[...]
    ["evento_terminale"]``, vedi nota di modulo.
    """
    risultati = []
    riepilogo = {}

    # --- Stadio 1: riuso diretto di guida.integra_gravity_turn (Fase A
    # verticale + kick + Fase B gravity turn). m_vuoto e' la soglia di
    # burnout di QUESTO stadio, non la massa a vuoto finale del razzo.
    stadio_1 = stadi[0]
    m_vuoto_1 = m0 - stadio_1["m_prop"]

    risultato_stadio_1 = guida.integra_gravity_turn(
        m0,
        m_vuoto_1,
        stadio_1["spinta"],
        stadio_1["isp"],
        cd,
        area,
        v_kick,
        kick_angle_deg,
        h0=h0,
        t_max=t_max,
        v_min=v_min,
    )
    risultato_gt_1 = risultato_stadio_1["gravity_turn"]

    risultati.append(risultato_stadio_1)

    stato_effettivo = risultato_gt_1.y[:, -1]  # [x, h, v, gamma, m]
    massa_burnout_effettiva = stato_effettivo[4]
    nome_evento_1, t_evento_1 = _evento_terminale(risultato_gt_1)
    riepilogo["stadio_1"] = {
        "massa_burnout_effettiva": massa_burnout_effettiva,
        "massa_burnout_nominale": m_vuoto_1,
        "differenza_massa": massa_burnout_effettiva - m_vuoto_1,
        "t_burnout": t_evento_1,
        "evento_terminale": nome_evento_1,
    }

    # Massa di partenza dello stadio successivo: massa EFFETTIVA al
    # burnout meno la massa strutturale espulsa dello stadio appena
    # concluso (addendum reviewer, punto 3).
    massa_inizio_segmento = massa_burnout_effettiva - stadio_1["m_strut"]
    x_eff, h_eff, v_eff, gamma_eff = stato_effettivo[:4]

    # --- Stadi 2..N: solve_ivp diretto, riuso di
    # guida.derivate_stato_gravity_turn e dei suoi eventi difensivi.
    for indice_stadio_0based, stadio in enumerate(stadi[1:], start=1):
        numero_stadio = indice_stadio_0based + 1  # 1-based per leggibilita'

        m_vuoto_i = massa_inizio_segmento - stadio["m_prop"]
        mdot_i = stadio["spinta"] / (stadio["isp"] * G0)

        y0_i = np.array([x_eff, h_eff, v_eff, gamma_eff, massa_inizio_segmento])
        args_i = (stadio["spinta"], stadio["isp"], cd, area, mdot_i)

        evento_fine_prop_i = functools.partial(
            guida.evento_fine_propellente_2d, m_vuoto=m_vuoto_i
        )
        evento_fine_prop_i.terminal = True
        evento_fine_prop_i.direction = -1

        evento_vmin_i = functools.partial(guida.evento_velocita_minima, v_min=v_min)
        evento_vmin_i.terminal = True
        evento_vmin_i.direction = -1

        risultato_i = solve_ivp(
            guida.derivate_stato_gravity_turn,
            t_span=(0.0, t_max),
            y0=y0_i,
            method="RK45",
            args=args_i,
            events=(
                evento_fine_prop_i,
                guida.evento_impatto_suolo_2d,
                guida.evento_traiettoria_invalida,
                evento_vmin_i,
            ),
            rtol=1e-8,
            # Stesse scale/tolleranze di Fase B in guida.py (stesso
            # vettore di stato a 5 componenti [x, h, v, gamma, m]).
            atol=[1e-3, 1e-3, 1e-6, 1e-9, 1e-3],
            dense_output=False,
        )
        risultati.append(risultato_i)

        stato_effettivo = risultato_i.y[:, -1]
        massa_burnout_effettiva = stato_effettivo[4]
        nome_evento_i, t_evento_i = _evento_terminale(risultato_i)
        riepilogo[f"stadio_{numero_stadio}"] = {
            "massa_burnout_effettiva": massa_burnout_effettiva,
            "massa_burnout_nominale": m_vuoto_i,
            "differenza_massa": massa_burnout_effettiva - m_vuoto_i,
            "t_burnout": t_evento_i,
            "evento_terminale": nome_evento_i,
        }

        massa_inizio_segmento = massa_burnout_effettiva - stadio["m_strut"]
        x_eff, h_eff, v_eff, gamma_eff = stato_effettivo[:4]

    return {"risultati": risultati, "riepilogo": riepilogo}
