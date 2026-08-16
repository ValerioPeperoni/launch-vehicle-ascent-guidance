"""Guida attiva, fase atmosferica: gravity turn (Step 3).

Dopo un breve tratto di ascesa verticale pura (Fase A, fisicamente
identica allo Step 2 e riusata cosi' com'e' da ``lanciatore.dinamica``),
si applica un kick angle iniziale e si passa a un moto 2D nel piano
verticale (Fase B) in cui la spinta e' sempre allineata alla velocita'
(angolo di attacco nullo per ipotesi) e l'assetto (``gamma``, angolo di
rotta rispetto all'orizzonte locale — flight-path angle) evolve secondo
l'equazione classica del gravity turn.

Vettore di stato Fase B: ``y = [x, h, v, gamma, m]`` (distanza a terra,
quota, modulo velocita', flight-path angle in RADIANTI, massa). La
Fase A resta nello stato 1D ``[h, v, m]`` gia' definito in
``lanciatore.dinamica`` (nessuna riformulazione necessaria li').

Equazioni del moto (Fase B), spinta sempre allineata alla velocita'
(nessuna componente di portanza/forza normale: angolo di attacco nullo
per costruzione, coerente con l'ipotesi di gravity turn dichiarata in
CLAUDE.md)::

    dv/dt     = (T - D) / m - g(h) * sin(gamma)
    dgamma/dt = -(g(h) / v) * cos(gamma)
    dh/dt     = v * sin(gamma)
    dx/dt     = v * cos(gamma)
    dm/dt     = -mdot

con ``D = 0.5 * rho(h) * CD * area * v**2`` (drag lungo la velocita',
verso opposto; ``v > 0`` per costruzione in Fase B, la Fase B non parte
mai da velocita' nulla) e ``g(h)`` calcolata da
``lanciatore.gravita.accelerazione_gravita`` (riusata, non ridefinita).

Fonte esplicita dell'equazione di guida (obbligatoria, vincolo
CLAUDE.md "Guida, fase atmosferica ... riferimento Culler et al. 1957"):
Culler, G.J. e Fried, B.D., "General Theory of Gravity Turn
Trajectories", Journal of Applied Physics / ARS Journal, 1957. E' la
fonte dell'equazione ``dgamma/dt = -(g/v)*cos(gamma)``, derivata
imponendo che l'angolo di attacco sia nullo (spinta e velocita' sempre
allineate) durante tutta la fase atmosferica propulsa.

Nota sull'approssimazione geometrica (segnalata dal reviewer,
esplicitata qui come richiesto): ai fini del campo gravitazionale,
``h`` e' trattata come la distanza radiale dal centro Terra che entra
in ``g(h) = MU_TERRA/(R_TERRA+h)**2`` (gravita' sempre "verticale" nel
senso locale), mentre ``x`` resta una coordinata di piano piatto, non
curva. Questa e' l'approssimazione 2D non-rotante gia' dichiarata in
CLAUDE.md; nello Step 2 (moto puramente radiale) era irrilevante perche'
non esisteva una coordinata orizzontale, qui diventa un'assunzione
aggiuntiva esplicita: la curvatura terrestre nel legame tra ``x`` e la
geometria del campo di gravita' resta fuori scope (nessun accoppiamento
x-g introdotto).

Nota sulla non-sovrapposizione degli eventi difensivi di Fase B
(segnalata dal reviewer, nessuna modifica richiesta, resa esplicita qui
per un lettore futuro): ``evento_traiettoria_invalida`` (``gamma <= 0``)
ed ``evento_impatto_suolo_2d`` (``h <= 0``) non possono scattare "nello
stesso istante in conflitto" perche' ``dh/dt = v*sin(gamma)``: finche'
``gamma`` resta in ``(0, 90°]`` il seno e' positivo e la quota non puo'
scendere, quindi ``gamma`` deve necessariamente attraversare lo zero
(facendo scattare l'evento di traiettoria invalida) *prima* che ``h``
possa iniziare a diminuire verso l'evento di impatto suolo.
"""

import functools

import numpy as np
from scipy.integrate import solve_ivp

from lanciatore import dinamica
from lanciatore.atmosfera import densita
from lanciatore.costanti import G0
from lanciatore.gravita import accelerazione_gravita


def derivate_stato_gravity_turn(t, y, spinta, isp, cd, area, mdot):
    """Derivate del vettore di stato ``[x, h, v, gamma, m]`` (Fase B, gravity turn).

    Equazioni (vedi docstring di modulo per fonte e note):
        dv/dt     = (T - D) / m - g(h) * sin(gamma)
        dgamma/dt = -(g(h) / v) * cos(gamma)
        dh/dt     = v * sin(gamma)
        dx/dt     = v * cos(gamma)
        dm/dt     = -mdot

    Riusa senza reimplementarle ``lanciatore.atmosfera.densita`` e
    ``lanciatore.gravita.accelerazione_gravita``.

    Parametri
    ---------
    t : float
        Tempo (non usato esplicitamente nelle equazioni, richiesto
        dalla firma di solve_ivp).
    y : array_like, shape (5,)
        Stato corrente ``[x, h, v, gamma, m]``, con ``gamma`` in
        radianti.
    spinta : float
        Spinta del motore, N. Sempre allineata alla velocita' (angolo
        di attacco nullo per ipotesi di gravity turn).
    isp : float
        Impulso specifico, s. Non usato direttamente qui (mdot e'
        gia' calcolata a monte con G0 costante), presente in firma
        solo per coerenza con lo stile gia' adottato in
        ``lanciatore.dinamica.derivate_stato``.
    cd : float
        Coefficiente di resistenza aerodinamica, costante.
    area : float
        Area di riferimento (sezione trasversale), m^2.
    mdot : float
        Portata massica, kg/s (positiva; la massa diminuisce).

    Ritorna
    -------
    numpy.ndarray, shape (5,)
        ``[dx/dt, dh/dt, dv/dt, dgamma/dt, dm/dt]``
    """
    x, h, v, gamma, m = y

    g = accelerazione_gravita(h)
    rho = densita(h)
    drag = 0.5 * rho * cd * area * v**2

    dv_dt = (spinta - drag) / m - g * np.sin(gamma)
    dgamma_dt = -(g / v) * np.cos(gamma)
    dh_dt = v * np.sin(gamma)
    dx_dt = v * np.cos(gamma)
    dm_dt = -mdot

    return np.array([dx_dt, dh_dt, dv_dt, dgamma_dt, dm_dt])


# ---------------------------------------------------------------------
# Eventi Fase B (stato a 5 componenti [x, h, v, gamma, m]).
#
# Reimplementati da zero rispetto agli eventi 1D dello Step 2 perche'
# l'indice della massa/quota nel vettore di stato e' diverso (5
# componenti vs 3): NON e' duplicazione della fisica di drag/gravita'/
# massa (che resta unica, in atmosfera.py/gravita.py/nelle derivate
# sopra), e' la replicazione minima e necessaria del pattern "funzione
# evento con terminal/direction come attributi" per un vettore di stato
# diverso, gia' stabilito nello Step 2.
# ---------------------------------------------------------------------


def evento_impatto_suolo_2d(t, y, spinta, isp, cd, area, mdot):
    """Evento di sicurezza Fase B: si annulla se la quota torna a zero (o sotto).

    Analogo a ``lanciatore.dinamica.evento_impatto_suolo`` (Step 2), qui
    per lo stato a 5 componenti. Vedi la docstring di modulo per la nota
    sulla non-sovrapposizione con ``evento_traiettoria_invalida``.
    """
    return y[1]


evento_impatto_suolo_2d.terminal = True
evento_impatto_suolo_2d.direction = -1


def evento_traiettoria_invalida(t, y, spinta, isp, cd, area, mdot):
    """Evento di sicurezza Fase B: si attiva se ``gamma`` scende a 0 o sotto.

    L'equilibrio del gravity turn puro a ``gamma = 90 gradi`` e'
    instabile (vedi test dedicato in tests/test_guida.py): un kick
    angle troppo grande o parametri veicolo poco realistici possono far
    collassare ``gamma`` verso il basso fino a una traiettoria che
    punta sempre piu' verso l'orizzonte ("tuffo"). Questo evento
    intercetta il caso e lo rende un fallimento esplicito e
    diagnosticabile invece di un'integrazione silenziosamente
    patologica.
    """
    return y[3]


evento_traiettoria_invalida.terminal = True
evento_traiettoria_invalida.direction = -1


def evento_fine_propellente_2d(t, y, spinta, isp, cd, area, mdot, m_vuoto):
    """Evento terminale Fase B: si annulla quando la massa raggiunge m_vuoto.

    Analogo a ``lanciatore.dinamica.evento_fine_propellente`` (Step 2).

    Nota implementativa: ``derivate_stato_gravity_turn`` (equazioni del
    punto 1.2 del piano) NON prende ``m_vuoto`` come parametro (non
    serve alla fisica del moto), quindi ``m_vuoto`` non fa parte della
    tupla ``args`` condivisa con la funzione derivate passata a
    ``solve_ivp``. Va percio' legato con ``functools.partial`` PRIMA di
    passare questa funzione a ``solve_ivp`` (stessa tecnica gia'
    richiesta per ``evento_velocita_kick``/v_kick in Fase A e per
    ``evento_velocita_minima``/v_min qui sotto — vedi
    ``integra_gravity_turn``), con ``terminal``/``direction`` impostati
    sull'oggetto ``functools.partial`` effettivamente passato, non su
    questa funzione "template".
    """
    return y[4] - m_vuoto


def evento_velocita_minima(t, y, spinta, isp, cd, area, mdot, v_min):
    """Evento difensivo Fase B: si annulla se v scende sotto una soglia minima.

    Safeguard per la singolarita' ``1/v`` in ``dgamma/dt = -(g(h)/v) *
    cos(gamma)`` (addendum reviewer, punto 2): a differenza dello
    Step 2, qui compare un termine non protetto in ``1/v``. Se ``v``
    collassasse verso zero (per drag eccessivo o parametri non fisici)
    l'integrazione diventerebbe numericamente patologica prima che
    ``evento_traiettoria_invalida`` possa intercettarla in modo
    affidabile. Questo evento ferma l'integrazione in modo esplicito e
    diagnosticabile, analogo nello spirito a
    ``evento_traiettoria_invalida`` per ``gamma``.

    Stessa nota implementativa di ``evento_fine_propellente_2d``: v_min
    va legato con ``functools.partial`` prima di passare la funzione a
    ``solve_ivp``, con ``terminal``/``direction`` impostati sul
    ``functools.partial`` risultante.
    """
    return y[2] - v_min


def evento_velocita_kick(t, y, spinta, isp, cd, area, mdot, m_vuoto, v_kick):
    """Evento terminale Fase A: si attiva quando la velocita' raggiunge v_kick.

    Fase A usa ``lanciatore.dinamica.derivate_stato`` cosi' com'e'
    (stato ``[h, v, m]``, firma ``args = (spinta, isp, cd, area, mdot,
    m_vuoto)`` gia' fissata dallo Step 2): ``v_kick`` non fa parte di
    quella tupla condivisa, quindi va legato con ``functools.partial``
    prima di passare questa funzione a ``solve_ivp`` (vedi
    ``integra_gravity_turn``), con ``terminal=True`` e ``direction=1``
    impostati sull'oggetto risultante dal partial (la velocita' cresce
    da 0 verso v_kick durante la Fase A, quindi l'attraversamento e' in
    direzione crescente).
    """
    return y[1] - v_kick


def integra_gravity_turn(
    m0,
    m_vuoto,
    spinta,
    isp,
    cd,
    area,
    v_kick,
    kick_angle_deg,
    h0=0.0,
    t_max=1000.0,
    v_min=1.0,
):
    """Integra ascesa verticale (Fase A) + gravity turn (Fase B).

    Fase A: ascesa verticale pura, identica per costruzione allo
    Step 2, riusando direttamente ``lanciatore.dinamica.derivate_stato``
    e gli eventi difensivi gia' definiti li' (``evento_fine_propellente``,
    ``evento_impatto_suolo``), insieme al nuovo evento locale
    ``evento_velocita_kick``. Integra fino a quando ``v`` raggiunge
    ``v_kick`` (oppure, in caso patologico, fino a uno degli altri due
    eventi difensivi, nel qual caso viene sollevato un errore esplicito
    — vedi sotto).

    Kick (istantaneo, a fine Fase A): lo stato di Fase B riparte da
    ``x0 = 0``, ``h0 = h_kick``, ``v0 = v_kick_effettivo``,
    ``gamma0 = radians(90 - kick_angle_deg)``, ``m0 = m_kick``, dove
    ``h_kick``, ``v_kick_effettivo`` e ``m_kick`` sono presi dallo
    STATO EFFETTIVO dell'evento di Fase A (``risultato_a.y[:, -1]``),
    non dal parametro nominale ``v_kick`` ne' da un ricalcolo separato:
    e' la stessa fonte (stato dell'evento) gia' usata per h0/m0, per
    evitare un disallineamento con la tolleranza del root-finder di
    ``solve_ivp`` (addendum reviewer, punto 4).

    Fase B: gravity turn vero e proprio (``derivate_stato_gravity_turn``)
    fino all'evento "fine propellente" dello stadio corrente, oppure
    fino a uno dei tre eventi difensivi (impatto suolo, traiettoria
    invalida, velocita' minima).

    Parametri
    ---------
    m0, m_vuoto, spinta, isp, cd, area : float
        Parametri veicolo, stesso significato dello Step 2
        (``lanciatore.dinamica.integra_ascesa_verticale``).
    v_kick : float
        Velocita' verticale (Fase A) a cui applicare il kick angle,
        m/s. Parametro di design della legge di guida, non una
        costante fisica di progetto (resta un parametro esplicito
        della funzione, nessun default hard-coded).
    kick_angle_deg : float
        Kick angle, in GRADI (convertito internamente in radianti):
        ``gamma0 = 90 - kick_angle_deg``. Parametro di design della
        legge di guida.
    h0 : float, opzionale
        Quota iniziale di Fase A, m (default 0.0).
    t_max : float, opzionale
        Tempo massimo di integrazione per ciascuna fase, s (default
        1000.0): limite di sicurezza, entrambe le fasi si fermano
        prima per uno degli eventi terminali in condizioni normali.
    v_min : float, opzionale
        Soglia di velocita' minima per l'evento difensivo
        ``evento_velocita_minima`` in Fase B, m/s (default 1.0):
        piccola a sufficienza da non interferire con un gravity turn
        nominale (dove v cresce durante tutta la fase propulsa), ma
        non nulla, per fermare l'integrazione prima che ``1/v`` in
        ``dgamma/dt`` diventi numericamente patologico.

    Ritorna
    -------
    dict
        ``{"fase_verticale": risultato_a, "gravity_turn": risultato_b}``,
        entrambi oggetti risultato grezzi di ``solve_ivp`` (nessun
        risultato nascosto o post-processato: il chiamante/test deve
        poter accedere direttamente a ``.t``, ``.y``, ``.t_events``,
        ``.status`` di ciascuna fase).

    Solleva
    -------
    ValueError
        Se ``spinta <= m0 * g(h0)``: il lift-off non decolla per
        definizione fisica. Stesso controllo, stessa formula e stesso
        stile di messaggio di
        ``lanciatore.dinamica.integra_ascesa_verticale`` (Step 2).
        Nota: il controllo e' qui necessariamente RIPETUTO (non
        importato) perche' in ``dinamica.py`` e' inline dentro
        ``integra_ascesa_verticale`` e non una funzione standalone
        importabile, e questo step non modifica ``dinamica.py`` (vedi
        piano, punto 1.6) — la ripetizione e' quindi esplicita e
        commentata come mirror intenzionale del controllo di Step 2,
        non una duplicazione silenziosa della logica.
    RuntimeError
        Se l'evento "raggiunta velocita' di kick" non si attiva prima
        di ``t_max`` o prima di uno degli altri due eventi difensivi di
        Fase A (fine propellente/impatto suolo): configurazione non
        valida per questo step, da segnalare esplicitamente (addendum
        reviewer, punto 1c).
    """
    # --- Controllo di liftoff (mirror esplicito del controllo in
    # lanciatore.dinamica.integra_ascesa_verticale, vedi nota "Solleva"
    # sopra: dinamica.py non va modificato in questo step, quindi il
    # controllo non e' importabile come funzione standalone e viene
    # ripetuto qui, consapevolmente, non duplicato in modo silenzioso).
    g0_quota_iniziale = accelerazione_gravita(h0)
    if spinta <= m0 * g0_quota_iniziale:
        raise ValueError(
            "Spinta insufficiente al lift-off: spinta="
            f"{spinta!r} N <= m0*g(h0)={m0 * g0_quota_iniziale!r} N. "
            "Il veicolo non decolla con questi parametri."
        )

    # Portata massica: definita con G0 costante (convenzione
    # internazionale per Isp), condivisa da entrambe le fasi (il
    # bruciamento e' continuo attraverso il kick, che e' istantaneo).
    mdot = spinta / (isp * G0)

    # --- Fase A: ascesa verticale pura, riuso diretto di dinamica.py.
    y0_a = np.array([h0, 0.0, m0])
    args_a = (spinta, isp, cd, area, mdot, m_vuoto)

    evento_kick = functools.partial(evento_velocita_kick, v_kick=v_kick)
    evento_kick.terminal = True
    evento_kick.direction = 1

    risultato_a = solve_ivp(
        dinamica.derivate_stato,
        t_span=(0.0, t_max),
        y0=y0_a,
        method="RK45",
        args=args_a,
        # Eventi difensivi di dinamica.py inclusi ANCHE qui (addendum
        # reviewer, punto 1b), insieme al nuovo evento di kick: la Fase
        # A non deve bypassare silenziosamente le protezioni gia'
        # stabilite nello Step 2.
        events=(
            evento_kick,
            dinamica.evento_fine_propellente,
            dinamica.evento_impatto_suolo,
        ),
        rtol=1e-8,
        atol=[1e-3, 1e-6, 1e-3],
        dense_output=False,
    )

    t_events_kick = risultato_a.t_events[0]
    if len(t_events_kick) == 0:
        raise RuntimeError(
            "Velocita' di kick (v_kick="
            f"{v_kick!r} m/s) non raggiunta entro t_max={t_max!r} s "
            "ne' prima degli eventi difensivi di Fase A (fine "
            "propellente/impatto suolo): configurazione non valida "
            "per il gravity turn con questi parametri."
        )

    # Stato EFFETTIVO all'evento di kick (addendum reviewer, punto 4):
    # h_kick, v_kick_effettivo e m_kick sono presi dall'ultimo punto
    # della soluzione di Fase A (che coincide con l'istante
    # dell'evento terminale), non dal parametro nominale v_kick.
    h_kick, v_kick_effettivo, m_kick = risultato_a.y[:, -1]

    # --- Kick: discontinuita' istantanea solo su gamma.
    gamma0 = np.radians(90.0 - kick_angle_deg)
    y0_b = np.array([0.0, h_kick, v_kick_effettivo, gamma0, m_kick])

    # --- Fase B: gravity turn.
    args_b = (spinta, isp, cd, area, mdot)

    evento_fine_prop_2d = functools.partial(
        evento_fine_propellente_2d, m_vuoto=m_vuoto
    )
    evento_fine_prop_2d.terminal = True
    evento_fine_prop_2d.direction = -1

    evento_vmin = functools.partial(evento_velocita_minima, v_min=v_min)
    evento_vmin.terminal = True
    evento_vmin.direction = -1

    risultato_b = solve_ivp(
        derivate_stato_gravity_turn,
        t_span=(0.0, t_max),
        y0=y0_b,
        method="RK45",
        args=args_b,
        events=(
            evento_fine_prop_2d,
            evento_impatto_suolo_2d,
            evento_traiettoria_invalida,
            evento_vmin,
        ),
        rtol=1e-8,
        # Scale molto diverse tra le componenti di stato: x, h in
        # metri (~10^2-10^5), v in m/s (~10-10^3), gamma in radianti
        # (~10^-2-1, richiede atol molto piu' stretto), m in kg
        # (~10^3-10^5).
        atol=[1e-3, 1e-3, 1e-6, 1e-9, 1e-3],
        dense_output=False,
    )

    return {"fase_verticale": risultato_a, "gravity_turn": risultato_b}
