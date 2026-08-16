"""Dinamica a punto materiale, ascesa verticale pura, singolo stadio (Step 2).

Caso di test piu' semplice possibile prima di introdurre la guida
(Step 3/4): moto puramente verticale (1 grado di liberta' spaziale,
la quota h), un solo stadio, nessuna guida — l'angolo spinta-velocita'
e' fisso a 90 gradi rispetto all'orizzonte per costruzione.

Vettore di stato: y = [h, v, m] (quota, velocita' verticale, massa
corrente).

Equazioni del moto:
    dh/dt = v
    dv/dt = (T - D - m*g(h)) / m
    dm/dt = -mdot

con:
    D = 0.5 * rho(h) * CD * area * v^2 * sign(v)
        (drag sempre opposto al moto; per v >= 0 si riduce a
        D = 0.5*rho(h)*CD*area*v^2 con verso negativo. La forma con
        sign(v) e' scritta in modo generale per riusabilita' futura,
        quando dallo Step 3 in poi la velocita' potra' avere anche
        una componente non puramente verticale)
    mdot = T / (Isp * G0), costante durante il bruciamento.

Nota importante (da non confondere): la portata massica mdot usa
G0 (costante, definizione internazionale di Isp), mentre dv/dt usa
g(h) (variabile con la quota, vedi lanciatore/gravita.py). Sono due
grandezze fisicamente distinte anche se numericamente vicine a h=0.
"""

import numpy as np
from scipy.integrate import solve_ivp

from lanciatore.atmosfera import densita
from lanciatore.costanti import CD, G0
from lanciatore.gravita import accelerazione_gravita


def derivate_stato(t, y, spinta, isp, cd, area, mdot, m_vuoto):
    """Derivate del vettore di stato [h, v, m] per l'ascesa verticale pura.

    Parametri
    ---------
    t : float
        Tempo (non usato esplicitamente nelle equazioni, richiesto
        dalla firma di solve_ivp).
    y : array_like, shape (3,)
        Stato corrente [h, v, m].
    spinta : float
        Spinta del motore, N. Costante durante il bruciamento.
    isp : float
        Impulso specifico, s.
    cd : float
        Coefficiente di resistenza aerodinamica, costante.
    area : float
        Area di riferimento (sezione trasversale), m^2.
    mdot : float
        Portata massica, kg/s (positiva; la massa diminuisce). Calcolata
        una sola volta a monte con G0 costante (spinta/(isp*G0)), non
        ricalcolata qui ad ogni chiamata.
    m_vuoto : float
        Massa a vuoto, kg. Dopo l'evento di fine propellente la massa
        non deve scendere sotto questo valore; la funzione derivate
        continua comunque a restituire -mdot fino a quando l'evento
        terminale non ferma l'integrazione (vedi evento_fine_propellente).

    Ritorna
    -------
    numpy.ndarray, shape (3,)
        [dh/dt, dv/dt, dm/dt]
    """
    h, v, m = y

    g = accelerazione_gravita(h)
    rho = densita(h)
    drag = 0.5 * rho * cd * area * v**2 * np.sign(v)

    dh_dt = v
    dv_dt = (spinta - drag - m * g) / m
    dm_dt = -mdot

    return np.array([dh_dt, dv_dt, dm_dt])


def evento_fine_propellente(t, y, spinta, isp, cd, area, mdot, m_vuoto):
    """Evento terminale: si annulla quando la massa corrente raggiunge m_vuoto.

    Stessa firma delle derivate di stato (richiesto da solve_ivp quando
    si passano args comuni a fun ed eventi). direction=-1 e terminal=True
    sono impostati come attributi della funzione, non come parametri: e'
    il pattern richiesto da solve_ivp per gli eventi.

    Questo e' il pattern generale che verra' riusato per gli eventi di
    staging multiplo allo Step 5 (qui in forma "degenere" con un solo
    stadio): cambia solo il numero di eventi/stadi, non la struttura.
    """
    return y[2] - m_vuoto


evento_fine_propellente.terminal = True
evento_fine_propellente.direction = -1


def evento_impatto_suolo(t, y, spinta, isp, cd, area, mdot, m_vuoto):
    """Evento di sicurezza: si annulla se la quota torna a zero (o sotto).

    Con i parametri di test attesi (T/W0 > 1) questo evento non dovrebbe
    mai attivarsi durante un'ascesa propulsa valida: e' incluso come
    protezione difensiva contro parametri non fisici (es. spinta
    insufficiente passata per errore) e per stabilire fin da ora il
    pattern "piu' eventi terminali coesistono", necessario allo
    Step 5 con lo staging multiplo.

    Caso limite a t=0: con h0=0.0 il valore della funzione evento a
    t=0 e' 0, ma solve_ivp rileva un evento sul cambio di segno
    *durante* l'integrazione, non sul valore isolato all'istante
    iniziale. Con v0=0 e dv/dt(0) > 0 (spinta netta positiva per il
    controllo T > m0*g(0) fatto a monte), la quota si allontana da
    zero immediatamente dopo t=0 e l'evento non scatta spuriamente.
    """
    return y[0]


evento_impatto_suolo.terminal = True
evento_impatto_suolo.direction = -1


def integra_ascesa_verticale(
    m0, m_vuoto, spinta, isp, cd, area, h0=0.0, v0=0.0, t_max=1000.0
):
    """Integra l'ascesa verticale pura (singolo stadio, nessuna guida).

    Parametri
    ---------
    m0 : float
        Massa iniziale (con propellente), kg.
    m_vuoto : float
        Massa a vuoto (senza propellente), kg.
    spinta : float
        Spinta del motore, N. Costante.
    isp : float
        Impulso specifico, s.
    cd : float
        Coefficiente di resistenza aerodinamica, costante.
    area : float
        Area di riferimento (sezione trasversale), m^2.
    h0 : float, opzionale
        Quota iniziale, m (default 0.0, livello del mare).
    v0 : float, opzionale
        Velocita' verticale iniziale, m/s (default 0.0).
    t_max : float, opzionale
        Tempo massimo di integrazione, s (default 1000.0): limite di
        sicurezza, l'integrazione si ferma prima per uno degli eventi
        terminali (fine propellente o, in caso patologico, impatto
        suolo) in condizioni normali.

    Ritorna
    -------
    scipy.integrate.OdeResult
        Oggetto risultato grezzo di solve_ivp (non un sottoinsieme
        estratto): il chiamante/test deve poter accedere direttamente a
        result.t, result.y, result.t_events, result.status per
        verificare quale evento ha fermato l'integrazione.

    Solleva
    -------
    ValueError
        Se spinta <= m0 * g(0): il "lift-off" non decolla per
        definizione fisica (spinta netta iniziale non positiva), quindi
        l'integrazione non avrebbe senso fisico da compiere.
    """
    g0_quota_iniziale = accelerazione_gravita(h0)
    if spinta <= m0 * g0_quota_iniziale:
        raise ValueError(
            "Spinta insufficiente al lift-off: spinta="
            f"{spinta!r} N <= m0*g(h0)={m0 * g0_quota_iniziale!r} N. "
            "Il veicolo non decolla con questi parametri."
        )

    # Portata massica: definita con G0 costante (convenzione
    # internazionale per Isp), anche se la gravita' del moto usa g(h)
    # variabile. Calcolata una sola volta, non ad ogni chiamata delle
    # derivate.
    mdot = spinta / (isp * G0)

    y0 = np.array([h0, v0, m0])
    args = (spinta, isp, cd, area, mdot, m_vuoto)

    risultato = solve_ivp(
        derivate_stato,
        t_span=(0.0, t_max),
        y0=y0,
        method="RK45",
        args=args,
        events=(evento_fine_propellente, evento_impatto_suolo),
        rtol=1e-8,
        atol=[1e-3, 1e-6, 1e-3],
        dense_output=False,
    )

    return risultato
