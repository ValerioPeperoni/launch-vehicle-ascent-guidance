"""Guida attiva, fase esoatmosferica: tangente lineare (Step 4).

Vettore di stato ``y = [x, h, vx, vh, m]`` (cartesiano, non polare come
la Fase B del gravity turn dello Step 3): qui il controllo ``theta(t)``
e' gia' una funzione esplicita ed elementare di ``t``, quindi non serve
ricostruire l'angolo della velocita' per applicare la legge di guida.

Equazioni del moto (fase esoatmosferica: nessun drag, nessuna densita'
atmosferica; gravita' linearizzata a valore costante sull'arco di
manovra -- vedi sotto)::

    dx/dt  = vx
    dh/dt  = vh
    dvx/dt = (spinta/m) * cos(theta(t))
    dvh/dt = (spinta/m) * sin(theta(t)) - g_costante
    dm/dt  = -mdot

con ``theta(t) = arctan(A + B*t)`` (legge di guida a tangente lineare) e
``g_costante`` un parametro passato esplicitamente dal chiamante (es.
valutato una volta con ``lanciatore.gravita.accelerazione_gravita`` alla
quota iniziale della fase), NON ricalcolato a ogni passo
dell'integrazione: la costanza di ``g`` e' l'approssimazione dichiarata
intrinseca alla tecnica della tangente lineare (e' cio' che la rende
"lineare" invece che accoppiata a g(h) variabile), non una scorciatoia
aggiuntiva di questo modulo.

Fonte della tecnica (vincolo CLAUDE.md, "Guida, fase esoatmosferica"):
Perkins, F.M., "Derivation of Linear-Tangent Steering Laws". La
derivazione qui sotto e' rifatta da principi primi (calcolo delle
variazioni / principio del massimo di Pontryagin), non riprodotta
letteralmente dal testo originale (nessun accesso web disponibile,
stesso limite gia' gestito per Culler & Fried allo Step 3).

Derivazione (principio del massimo di Pontryagin)
--------------------------------------------------
Problema di controllo ottimo: minimizzare il tempo per raggiungere una
velocita' terminale ``(vx_f, vh_f)`` assegnata, con modulo di spinta
costante, nessun drag (fase esoatmosferica), gravita' linearizzata a
valore costante sull'arco di manovra. Stato: posizione ``(x, h)`` e
velocita' ``(vx, vh)``. Hamiltoniano::

    H = 1 + lambda_x*vx + lambda_h*vh
          + lambda_vx*a(t)*cos(theta) + lambda_vh*(a(t)*sin(theta) - g)

con ``a(t) = spinta/m(t)``.

Equazioni dei costati: ``dlambda_x/dt = dlambda_h/dt = 0`` (costanti);
``dlambda_vx/dt = -lambda_x`` (costante) -> ``lambda_vx(t)`` lineare in
``t``; ``dlambda_vh/dt = -lambda_h`` (costante) -> ``lambda_vh(t)``
lineare in ``t``. Il controllo ottimo massimizza
``lambda_vx*cos(theta) + lambda_vh*sin(theta)``, quindi
``(cos(theta), sin(theta))`` parallelo a ``(lambda_vx(t), lambda_vh(t))``:
``tan(theta(t)) = lambda_vh(t) / lambda_vx(t)``.

Condizione di trasversalita' e ``lambda_x = 0``: la posizione
orizzontale finale ``x_f`` e' LIBERA (nessun vincolo su dove il
veicolo si trovi orizzontalmente all'istante finale, solo sulla quota e
sulla velocita'). La condizione di trasversalita' del principio del
massimo impone che il costato coniugato a una coordinata di stato finale
libera si annulli: ``lambda_x = 0`` segue quindi DIRETTAMENTE da questa
condizione al contorno, non e' una "scelta di sistema di riferimento"
arbitraria. Con ``lambda_x = 0``, ``dlambda_vx/dt = -lambda_x = 0``,
quindi ``lambda_vx(t)`` resta COSTANTE, mentre ``lambda_vh(t)`` resta
lineare in ``t`` (perche' ``lambda_h`` in generale NON e' nullo, vedi
sotto). Si ottiene percio' esattamente::

    tan(theta(t)) = A + B*t

con ``A = lambda_vh(0) / lambda_vx`` e ``B = -lambda_h / lambda_vx``, due
costanti libere determinate dalle condizioni al contorno sulla velocita'
finale (problema di determinazione di A, B per un target reale: FUORI
SCOPE di questo step, vedi sotto).

Perche' serve h_f vincolata (B != 0): per lo stesso ragionamento di
trasversalita' applicato alla quota, se ANCHE ``h_f`` fosse libera si
avrebbe ``lambda_h = 0`` identicamente (nessun vincolo da soddisfare su
``h`` all'istante finale), e quindi ``B = -lambda_h/lambda_vx = 0``: la
legge di guida collasserebbe a un angolo costante ``theta(t) = arctan(A)``
(niente "tangente lineare" vera). La tangente lineare vera (con ``B``
strutturalmente diverso da zero) richiede quindi ESPLICITAMENTE che la
quota finale ``h_f`` sia VINCOLATA (tipicamente il target di iniezione
orbitale, es. l'apogeo/quota di inserimento), mentre la distanza a terra
``x_f`` resta libera. Questa asimmetria (``x_f`` libera, ``h_f``
fissata) e' l'assunzione che rende l'intera derivazione autoconsistente
ed e' scritta qui esplicitamente, non lasciata implicita.

Convenzione dell'angolo theta
------------------------------
``theta`` e' l'angolo della SPINTA rispetto all'orizzontale locale, in
un riferimento CARTESIANO FISSO (non relativo alla velocita' del
veicolo). Questo e' concettualmente diverso da ``gamma`` dello Step 3
(gravity turn): li' ``gamma`` era l'angolo della VELOCITA' rispetto
all'orizzonte, con la spinta vincolata ad essere allineata alla
velocita' stessa (angolo di attacco nullo per ipotesi). Qui, al
contrario, la spinta puo' NON essere allineata alla velocita' -- e'
guida attiva vera (il controllo theta(t) e' una legge di guida
indipendente, non un vincolo di assetto passivo).

Nota sul dominio di arctan (nessuna modifica richiesta, da tenere
presente): ``arctan(A + B*t)`` restituisce sempre valori nell'intervallo
aperto ``(-90 gradi, 90 gradi)``, quindi ``cos(theta) > 0`` sempre: la
parametrizzazione ``tan(theta) = A + B*t`` impone IMPLICITAMENTE che la
componente orizzontale della spinta sia sempre non negativa (mai spinta
"all'indietro", nessuna inversione della componente orizzontale). E'
ragionevole per un'iniezione orbitale (si vuole sempre accelerare in
avanti orizzontalmente), ma e' un vincolo implicito della
parametrizzazione stessa, non un'ipotesi fisica aggiuntiva imposta a
mano -- va tenuto presente se in step futuri servisse un profilo con
spinta che punta all'indietro (fuori scope qui).

Confine di questo step (scope deciso esplicitamente)
------------------------------------------------------
Questo modulo implementa e verifica la LEGGE di guida (la dinamica sotto
``tan(theta) = A + B*t`` con ``A``, ``B`` DATI/noti), non la soluzione
del problema al contorno generale (trovare ``A``, ``B`` per centrare un
target orbitale reale -- in letteratura PEG questo e' un problema
accoppiato risolto per iterazione/shooting, dato che con deplezione di
massa reale l'accelerazione ``a(t)`` non e' costante e l'integrale non
ha piu' forma chiusa). Trovare ``A``, ``B`` per un target di missione
reale resta esplicitamente FUORI SCOPE qui (rimandato a quando servira'
davvero, Step 5/6).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def derivate_stato_tangente_lineare(t, y, A, B, spinta, mdot, g_costante):
    """Derivate del vettore di stato ``[x, h, vx, vh, m]`` (tangente lineare).

    Equazioni (vedi docstring di modulo per la derivazione e la fonte):
        dx/dt  = vx
        dh/dt  = vh
        dvx/dt = (spinta/m) * cos(theta(t))
        dvh/dt = (spinta/m) * sin(theta(t)) - g_costante
        dm/dt  = -mdot

    con ``theta(t) = arctan(A + B*t)``.

    Parametri
    ---------
    t : float
        Tempo, s (usato esplicitamente qui, a differenza delle derivate
        della Fase B del gravity turn, Step 3: theta dipende da t).
    y : array_like, shape (5,)
        Stato corrente ``[x, h, vx, vh, m]``.
    A, B : float
        Coefficienti della legge di guida a tangente lineare,
        ``tan(theta(t)) = A + B*t`` (``B`` in 1/s). Determinati a monte
        dalle condizioni al contorno sulla velocita' finale (fuori
        scope in questo modulo, vedi docstring di modulo).
    spinta : float
        Spinta del motore, N.
    mdot : float
        Portata massica, kg/s (positiva; la massa diminuisce). Puo'
        essere 0.0 (caso limite usato per il confronto in forma chiusa,
        vedi tests/test_guida_esoatmosferica.py).
    g_costante : float
        Accelerazione di gravita', m/s^2, valutata UNA VOLTA a monte
        (es. con lanciatore.gravita.accelerazione_gravita alla quota
        iniziale della fase) e passata come costante per l'intera
        integrazione: la costanza di g e' l'approssimazione dichiarata
        della tecnica della tangente lineare, non ricalcolata qui.

    Ritorna
    -------
    numpy.ndarray, shape (5,)
        ``[dx/dt, dh/dt, dvx/dt, dvh/dt, dm/dt]``
    """
    x, h, vx, vh, m = y

    theta = np.arctan(A + B * t)
    a = spinta / m

    dx_dt = vx
    dh_dt = vh
    dvx_dt = a * np.cos(theta)
    dvh_dt = a * np.sin(theta) - g_costante
    dm_dt = -mdot

    return np.array([dx_dt, dh_dt, dvx_dt, dvh_dt, dm_dt])


def integra_tangente_lineare(
    m0, spinta, mdot, A, B, g_costante, tf, x0=0.0, h0=0.0, vx0=0.0, vh0=0.0
):
    """Integra la dinamica a tangente lineare su un orizzonte temporale fisso.

    Nessun evento terminale: ``tf`` e' fisso e dato (coerente col
    confine di scope di questo step, vedi docstring di modulo -- non e'
    ancora presente ne' staging ne' una condizione di arresto legata al
    raggiungimento di un target orbitale, che appartiene a step
    successivi).

    Parametri
    ---------
    m0 : float
        Massa iniziale, kg.
    spinta : float
        Spinta del motore, N.
    mdot : float
        Portata massica, kg/s (puo' essere 0.0).
    A, B : float
        Coefficienti della legge di guida, vedi
        ``derivate_stato_tangente_lineare``.
    g_costante : float
        Accelerazione di gravita' costante per l'integrazione, m/s^2.
    tf : float
        Istante finale dell'integrazione, s (``t_span = (0.0, tf)``).
    x0, h0, vx0, vh0 : float, opzionali
        Condizioni iniziali di posizione e velocita' (default 0.0).

    Ritorna
    -------
    scipy.integrate.OdeResult
        Risultato grezzo di ``solve_ivp`` (nessun post-processing
        nascosto, stesso principio gia' seguito negli step precedenti):
        il chiamante/test accede direttamente a ``.t``, ``.y``,
        ``.success``.
    """
    y0 = np.array([x0, h0, vx0, vh0, m0])
    args = (A, B, spinta, mdot, g_costante)

    risultato = solve_ivp(
        derivate_stato_tangente_lineare,
        t_span=(0.0, tf),
        y0=y0,
        method="RK45",
        args=args,
        rtol=1e-8,
        # Scale diverse tra le componenti di stato: x, h in metri
        # (~10^2-10^4), vx, vh in m/s (~10-10^3), m in kg (~10^2-10^3).
        atol=[1e-3, 1e-3, 1e-6, 1e-6, 1e-3],
        dense_output=False,
    )

    return risultato


# ---------------------------------------------------------------------------
# Problema inverso (Step 6): trovare A, B della guida a tangente lineare che
# realizzano una velocita' terminale (vx_target, vh_target) assegnata, con
# tf FISSO e deplezione di massa reale (mdot != 0 in generale, quindi
# a(t) = spinta/m(t) NON e' costante e la forma chiusa dello Step 4
# (asinh/sqrt) non si applica -- vedi "Confine di questo step" nella
# docstring di modulo sopra).
# ---------------------------------------------------------------------------

# Tolleranza di residuo ASSOLUTA (m/s) sulle componenti di velocita' per
# accettare una soluzione di fsolve come "convergente". Verificata sia sul
# flag ier di fsolve sia sul residuo effettivo (vedi motivazione nella
# docstring di risolvi_coefficienti_tangente_lineare): il flag ier==1 da
# solo NON e' una garanzia sufficiente, fsolve puo' dichiararsi convergente
# con un residuo che in pratica e' comunque troppo grande per i nostri
# scopi (o, piu' raramente, il contrario).
TOLLERANZA_RESIDUO_FSOLVE = 1e-6

# Tolleranza (usata sia in senso assoluto sia relativo, via numpy.allclose)
# per verificare che il guess primario e il guess perturbato convergano
# sulla STESSA soluzione (A, B). Questa e' la protezione empirica contro il
# fenomeno delle radici multiple/spurie descritto sotto: se le due
# soluzioni non coincidono entro questa tolleranza, il guess fornito non e'
# saldamente dentro un unico bacino di attrazione e il risultato non e'
# affidabile.
TOLLERANZA_CONFRONTO_SOLUZIONI = 1e-6


def _residuo_velocita_finale(
    coefficienti,
    m0,
    spinta,
    mdot,
    g_costante,
    tf,
    vx_target,
    vh_target,
    x0,
    h0,
    vx0,
    vh0,
):
    """Funzione residuo per fsolve: scarto tra velocita' terminale ottenuta e target.

    Integra la dinamica ESISTENTE (``integra_tangente_lineare``, Step 4,
    riusata senza alcuna modifica) con i coefficienti ``A, B`` correnti e
    restituisce ``[vx(tf) - vx_target, vh(tf) - vh_target]``. Nessuna
    reimplementazione della fisica: questa funzione e' puro glue code tra
    ``fsolve`` e l'integratore gia' verificato allo Step 4.
    """
    A, B = coefficienti
    risultato = integra_tangente_lineare(
        m0, spinta, mdot, A, B, g_costante, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0
    )
    vx_finale = risultato.y[2, -1]
    vh_finale = risultato.y[3, -1]
    return np.array([vx_finale - vx_target, vh_finale - vh_target])


def risolvi_coefficienti_tangente_lineare(
    m0,
    spinta,
    mdot,
    g_costante,
    tf,
    vx_target,
    vh_target,
    A0,
    B0,
    x0=0.0,
    h0=0.0,
    vx0=0.0,
    vh0=0.0,
):
    """Risolve A, B della guida a tangente lineare per una velocita' terminale target (problema inverso, Step 6).

    Trova, per shooting (``scipy.optimize.fsolve``), i coefficienti
    ``A, B`` tali che integrando ``derivate_stato_tangente_lineare`` (via
    ``integra_tangente_lineare``, entrambe RIUSATE senza modifiche
    dall'implementazione dello Step 4) su ``[0, tf]`` si ottenga
    ``vx(tf) ~= vx_target`` e ``vh(tf) ~= vh_target``. Con ``mdot != 0``
    reale l'accelerazione ``a(t) = spinta/m(t)`` non e' costante, quindi
    non esiste una forma chiusa (a differenza del caso limite ``mdot = 0``
    gia' verificato allo Step 4 con asinh/sqrt) e serve necessariamente
    root-finding numerico.

    Fenomeno delle radici multiple/spurie (scoperto empiricamente,
    documentato qui invece di essere nascosto -- vedi STATUS.md, Ciclo 6,
    addendum)
    ---------------------------------------------------------------------
    Il sistema ``(A, B) -> (vx(tf), vh(tf))`` NON e' iniettivo: per la
    stessa coppia di target esistono (almeno) due radici distinte di
    ``fsolve``, una "vera" (fisicamente intesa, tipicamente con ``B`` dello
    stesso segno della variazione attesa dell'angolo di guida) e una
    "spuria" con ``B`` all'incirca di segno OPPOSTO -- una sorta di
    soluzione speculare: su un arco temporale breve, un profilo
    ``theta(t)`` leggermente crescente e uno leggermente decrescente
    possono produrre una velocita' media integrata quasi identica.
    Verificato empiricamente (non solo sospettato): con un guess neutro
    (``B0 = 0``) o di segno opposto al vero, ``fsolve`` converge quasi
    sempre alla radice spuria con residuo a precisione di macchina --
    cioe' un successo "silenzioso" nel senso del criterio di arresto di
    ``fsolve``, ma un risultato fisicamente sbagliato (errori relativi
    osservati fino al 200% su ``B``, col segno invertito). Un tentativo di
    costruire un guess automatico risolvendo la forma chiusa dello Step 4
    (caso ``mdot=0``) NON risolve il problema: converge anch'esso in modo
    incoerente sulla radice giusta o sbagliata a seconda del caso, perche'
    il problema non e' "quanto e' vicino il guess" ma "in quale bacino di
    attrazione (segno di B) si parte" -- informazione che una formula
    chiusa approssimata non puo' fornire in generale senza gia' conoscere
    la risposta.

    Perche' A0, B0 sono OBBLIGATORI (nessun default, nessun calcolo
    automatico)
    ---------------------------------------------------------------------
    Di conseguenza diretta di quanto sopra, questa funzione NON calcola un
    guess iniziale automaticamente: il chiamante deve fornire ``A0, B0``
    gia' nel bacino di attrazione fisicamente corretto (tipicamente noto
    dal contesto di missione -- es. il segno atteso della variazione
    dell'angolo di guida -- oppure una soluzione vicina gia' nota, come
    avviene nei cicli iterativi reali di guida esplicita: il PEG usato
    dallo Space Shuttle non riparte mai da un default generico, riusa la
    soluzione del ciclo di guida precedente). Non c'e' NESSUNA garanzia di
    convergenza globale: questa funzione risolve un problema locale (dato
    un buon guess, trova la radice in quel bacino), non un problema di
    ottimizzazione globale.

    Safeguard di robustezza (rileva l'ambiguita', non la risolve)
    ---------------------------------------------------------------------
    Dopo la convergenza dal guess fornito ``(A0, B0)``, la funzione ripete
    la risoluzione da un SECONDO guess ottenuto perturbando il primo entro
    lo stesso bacino atteso (``A0*1.3, B0*1.3`` -- stesso segno, quindi se
    il bacino e' quello giusto entrambi i tentativi devono convergere alla
    STESSA radice) e verifica che le due soluzioni coincidano entro
    tolleranza stretta. Se non coincidono, la funzione solleva
    ``RuntimeError``: NON tenta di scegliere una delle due o di fare la
    media, perche' un disaccordo tra i due tentativi e' il sintomo che il
    guess fornito non e' saldamente dentro un bacino di attrazione stabile
    (situazione in cui qualunque soluzione restituita silenziosamente
    sarebbe inaffidabile).

    Parametri
    ---------
    m0 : float
        Massa iniziale, kg.
    spinta : float
        Spinta del motore, N.
    mdot : float
        Portata massica, kg/s (puo' essere 0.0, vedi test companion non
        circolare).
    g_costante : float
        Accelerazione di gravita' costante per l'integrazione, m/s^2
        (stesso significato di ``integra_tangente_lineare``).
    tf : float
        Istante finale FISSO, s (non un'incognita di questo problema,
        coerente con lo scope gia' stabilito allo Step 4).
    vx_target, vh_target : float
        Velocita' orizzontale/verticale desiderata all'istante ``tf``,
        m/s.
    A0, B0 : float
        Guess iniziale OBBLIGATORIO per ``A`` e ``B`` (nessun default:
        vedi motivazione sopra). Deve essere gia' nel bacino di
        attrazione fisicamente corretto.
    x0, h0, vx0, vh0 : float, opzionali
        Condizioni iniziali di posizione e velocita' (default 0.0), stesso
        significato di ``integra_tangente_lineare``.

    Ritorna
    -------
    (A, B) : tuple of float
        Coefficienti della legge di guida a tangente lineare che
        realizzano il target entro la tolleranza di residuo dichiarata.

    Solleva
    -------
    RuntimeError
        Se ``fsolve`` non converge (``ier != 1``) dal guess primario o dal
        guess perturbato, se il residuo effettivo (``info['fvec']``,
        non il solo flag ``ier``) supera ``TOLLERANZA_RESIDUO_FSOLVE`` per
        uno dei due tentativi, oppure se le due soluzioni (guess primario
        e guess perturbato) non coincidono entro
        ``TOLLERANZA_CONFRONTO_SOLUZIONI``.
    """

    def _risolvi_da_guess(A_guess, B_guess):
        guess = np.array([A_guess, B_guess], dtype=float)
        soluzione, info, ier, messaggio = fsolve(
            _residuo_velocita_finale,
            x0=guess,
            args=(m0, spinta, mdot, g_costante, tf, vx_target, vh_target, x0, h0, vx0, vh0),
            full_output=True,
        )
        residuo_max = float(np.max(np.abs(info["fvec"])))
        return soluzione, ier, residuo_max, messaggio

    soluzione_primaria, ier_primaria, residuo_primario, msg_primaria = _risolvi_da_guess(A0, B0)

    A_pert, B_pert = A0 * 1.3, B0 * 1.3
    soluzione_perturbata, ier_perturbata, residuo_perturbato, msg_perturbata = _risolvi_da_guess(
        A_pert, B_pert
    )

    if ier_primaria != 1 or residuo_primario > TOLLERANZA_RESIDUO_FSOLVE:
        raise RuntimeError(
            "risolvi_coefficienti_tangente_lineare: fsolve non converge dal "
            f"guess primario (A0={A0}, B0={B0}): ier={ier_primaria}, "
            f"residuo_max={residuo_primario:.3e} (tolleranza "
            f"{TOLLERANZA_RESIDUO_FSOLVE:.1e}). Messaggio scipy: {msg_primaria}"
        )

    if ier_perturbata != 1 or residuo_perturbato > TOLLERANZA_RESIDUO_FSOLVE:
        raise RuntimeError(
            "risolvi_coefficienti_tangente_lineare: fsolve non converge dal "
            f"guess perturbato (A0={A_pert}, B0={B_pert}): ier={ier_perturbata}, "
            f"residuo_max={residuo_perturbato:.3e} (tolleranza "
            f"{TOLLERANZA_RESIDUO_FSOLVE:.1e}). Messaggio scipy: {msg_perturbata}"
        )

    if not np.allclose(
        soluzione_primaria,
        soluzione_perturbata,
        rtol=TOLLERANZA_CONFRONTO_SOLUZIONI,
        atol=TOLLERANZA_CONFRONTO_SOLUZIONI,
    ):
        raise RuntimeError(
            "risolvi_coefficienti_tangente_lineare: il guess primario "
            f"(A0={A0}, B0={B0}) e il guess perturbato (A0={A_pert}, "
            f"B0={B_pert}) convergono su soluzioni diverse: primario "
            f"(A,B)={tuple(soluzione_primaria)}, perturbato "
            f"(A,B)={tuple(soluzione_perturbata)}. Questo e' il sintomo del "
            "fenomeno delle radici multiple/spurie documentato nella "
            "docstring di questa funzione: il guess fornito non e' "
            "saldamente dentro un unico bacino di attrazione, il risultato "
            "non e' affidabile."
        )

    A_soluzione, B_soluzione = soluzione_primaria
    return float(A_soluzione), float(B_soluzione)
