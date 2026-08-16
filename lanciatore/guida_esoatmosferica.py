"""Guida attiva, fase esoatmosferica: tangente lineare (Step 4/6, corretto Ciclo 7).

Vettore di stato ``y = [x, h, vx, vh, m]`` (cartesiano, non polare come
la Fase B del gravity turn dello Step 3): qui il controllo ``theta(t)``
e' gia' una funzione esplicita ed elementare di ``t``, quindi non serve
ricostruire l'angolo della velocita' per applicare la legge di guida.

Equazioni del moto (fase esoatmosferica: nessun drag, nessuna densita'
atmosferica)::

    dx/dt  = vx
    dh/dt  = vh
    dvx/dt = (spinta/m) * cos(theta(t)) - vh*vx/(R_TERRA+h)
    dvh/dt = (spinta/m) * sin(theta(t)) - g(h) + vx**2/(R_TERRA+h)
    dm/dt  = -mdot

con ``theta(t) = arctan(A + B*t)`` (legge di guida a tangente lineare) e
``g(h) = lanciatore.gravita.accelerazione_gravita(h)`` (riusata,
ricalcolata a ogni passo -- NON piu' una costante congelata a monte,
vedi sezione "Correzione fisica" sotto per il motivo).

Fonte della tecnica di guida (vincolo CLAUDE.md, "Guida, fase
esoatmosferica"): Perkins, F.M., "Derivation of Linear-Tangent Steering
Laws". La derivazione del controllo ottimo (sezione piu' sotto) e'
rifatta da principi primi (calcolo delle variazioni / principio del
massimo di Pontryagin), non riprodotta letteralmente dal testo
originale (nessun accesso web disponibile, stesso limite gia' gestito
per Culler & Fried allo Step 3).

Correzione fisica (Ciclo 7): il termine di sollievo centripeto
------------------------------------------------------------------
**La scoperta.** Assemblando per la prima volta una traiettoria
completa con dati reali (Falcon 9 Block 5, Wikipedia, vedi
``lanciatore/validazione.py`` e STATUS.md Ciclo 7), con il modello
originale a gravita' CONGELATA costante (``g_costante`` valutato una
volta a monte, mai piu' ricalcolato) la velocita' orbitale a 200 km
(~7788 m/s) risultava IRRAGGIUNGIBILE dallo Stadio 2 per qualunque
coppia ``A, B``, con un deficit di ~350 m/s anche a payload zero.
Causa: vicino alla velocita' orbitale il "sollievo centripeto" dovuto
alla curvatura della traiettoria attorno alla Terra, ``vx**2/(R_TERRA+h)``,
diventa comparabile alla gravita' stessa (per definizione, all'orbita
circolare ``v_orbitale**2/r = g(r)`` esattamente): un modello che lo
ignora sovrastima sistematicamente il "peso" da vincere man mano che
``vx`` cresce, e diventa via via piu' sbagliato proprio nel regime che
piu' conta per l'iniezione orbitale.

**La derivazione completa (coordinate polari locali).** Sia ``r =
R_TERRA + h`` la distanza dal centro della Terra e ``phi`` l'angolo
sotteso (downrange). Le equazioni del moto in campo centrale (nessuna
spinta, nessun drag) sono, in coordinate polari::

    r'' - r*(phi')**2 = -g(r)                    (radiale)
    r*phi'' + 2*r'*phi' = 0                        (tangenziale, <=>
                                                     d/dt[r**2 * phi'] = 0,
                                                     conservazione del
                                                     momento angolare)

Con ``vh = dh/dt = r'`` (velocita' radiale) e ``vx = r*phi'`` (velocita'
tangenziale, cio' che nel resto del modulo si chiama semplicemente "la
componente orizzontale della velocita'"):

- ``dvx/dt = d(r*phi')/dt = r'*phi' + r*phi'' = (vh/r)*vx + r*phi''``.
  Dall'equazione tangenziale, ``r*phi'' = -2*r'*phi' = -2*vh*vx/r``,
  quindi ``dvx/dt = vh*vx/r - 2*vh*vx/r = -vh*vx/r``.
- ``dvh/dt = r'' = r*(phi')**2 - g(r) = vx**2/r - g(r)``.

Aggiungendo la spinta (proiettata con l'angolo ``theta`` nel
riferimento cartesiano locale, vedi sezione "Convenzione dell'angolo
theta" sotto) si ottengono esattamente le due equazioni riportate in
cima a questo modulo, con ``r = R_TERRA + h``. **I due termini nuovi
(``-vh*vx/(R_TERRA+h)`` in ``dvx/dt`` e ``+vx**2/(R_TERRA+h)`` in
``dvh/dt``) nascono dalla STESSA trasformazione geometrica**: non sono
due correzioni indipendenti di cui si puo' scegliere di applicarne una
sola, sono i due lati della stessa identita'.

**La dimostrazione che serve ENTRAMBI i termini (non solo uno).** Nel
limite balistico (``spinta=0``), dalle equazioni sopra si verifica per
sostituzione diretta che sono conservati ESATTAMENTE:

- il momento angolare specifico ``L = r*vx = vx*(R_TERRA+h)``:
  ``dL/dt = r'*vx + r*dvx/dt = vh*vx + r*(-vh*vx/r) = vh*vx - vh*vx = 0``;
- l'energia meccanica specifica
  ``E = (vx**2+vh**2)/2 - MU_TERRA/(R_TERRA+h)``:
  ``dE/dt = vx*dvx/dt + vh*dvh/dt + (MU_TERRA/r**2)*vh``
  ``       = vx*(-vh*vx/r) + vh*(vx**2/r - g(r)) + (MU_TERRA/r**2)*vh``
  ``       = -vh*g(r) + vh*MU_TERRA/r**2 = 0`` (perche' ``g(r) =
  MU_TERRA/r**2`` esattamente, vedi ``lanciatore/gravita.py``).

Una prima proposta di correzione (Ciclo 7, bozza iniziale
dell'orchestratore, poi corretta dal reviewer) aggiungeva SOLO il
termine in ``dvh/dt`` (il "sollievo centripeto" piu' intuitivo,
motivato da ``v_orbitale**2/r = g``) senza il termine gemello in
``dvx/dt``. Verificato numericamente (``scipy.integrate.solve_ivp``,
``rtol=atol=1e-12``, condizioni iniziali arbitrarie, 200 s di
integrazione balistica): con la correzione COMPLETA (entrambi i
termini) ``L`` ed ``E`` si conservano a precisione di macchina (~1e-14
di variazione relativa); con la correzione a META' (solo ``dvh/dt``),
``E`` varia dello **0.197%** e ``L`` dell'**1.25%** sugli stessi 200 s
-- un errore fisico reale (violazione di conservazione), non un
dettaglio accademico. Questo modulo implementa la versione COMPLETA
(vedi test di conservazione dedicati in
``tests/test_guida_esoatmosferica.py``).

**Perche' ``g(h)`` ora si ricalcola a ogni passo, invece di restare
congelato.** Con il sollievo centripeto che dipende esplicitamente da
``h`` (tramite ``R_TERRA+h`` al denominatore) e la traiettoria che ora
sale di centinaia di km durante la fase esoatmosferica, congelare
anche ``g`` al valore iniziale reintrodurrebbe un errore sistematico
dello stesso tipo di quello appena corretto (anche se piu' piccolo:
``g(h)`` varia di pochi punti percentuali su poche centinaia di km,
contro il sollievo centripeto che vicino a ``v_orbitale`` e' comparabile
a ``g`` stesso) -- e non costa nulla in piu' (``accelerazione_gravita``
e' gia' vettorizzata e usata a ogni passo in Step 2/3/5). La firma di
``derivate_stato_tangente_lineare`` NON accetta piu' ``g_costante``:
e' un parametro rimosso, non un default nascosto.

**Asimmetria dichiarata col gravity turn (Step 3, NON toccato in
questo ciclo).** ``lanciatore/guida.py`` (Fase B del gravity turn)
continua a usare l'equazione originale di Culler & Fried,
``dgamma/dt = -(g(h)/v)*cos(gamma)``, senza alcun termine di sollievo
centripeto -- quella derivazione e' presa cosi' com'e' dalla fonte
citata allo Step 3 e non viene toccata da questa correzione. L'errore
introdotto da questa omissione e' stato quantificato con i dati reali
di validazione (Step 5/Ciclo 7, non stimato a occhio): a fine Stadio 1,
``v**2/(R_TERRA+h)`` vale il **12.2%** di ``g(h)`` con i dati di spinta
al livello del mare (v ~= 2756 m/s, h ~= 45.3 km) e il **17.3%** con i
dati in vuoto (v ~= 3266 m/s, h ~= 80 km stimata) -- non trascurabile,
ma sensibilmente piu' piccolo del quasi-100% raggiunto vicino alla
velocita' orbitale (dove il gap era stato scoperto in questo modulo).
La correzione al gravity turn resta esplicitamente FUORI SCOPE per
questo ciclo (nessuna modifica a ``lanciatore/guida.py``): l'asimmetria
e' dichiarata qui, non nascosta, e andra' rivalutata se in un ciclo
futuro il contributo di questa approssimazione al bilancio di delta-v
risultasse non trascurabile.

**Nota sulla derivazione del controllo ottimo (sezione sotto) rispetto
a questa correzione.** La derivazione a tangente lineare (Pontryagin,
sezione seguente) assume esplicitamente gravita' linearizzata a valore
costante sull'arco di manovra: e' quella linearizzazione a rendere il
problema risolvibile in forma chiusa con costati lineari in ``t``. I
due termini aggiunti qui rendono la dinamica EFFETTIVAMENTE integrata
piu' accurata di quella linearizzata usata per DERIVARE la legge di
controllo ``tan(theta) = A+B*t`` -- quindi, in senso stretto, con
questa correzione la legge a tangente lineare non e' piu' l'ottimo
matematico ESATTO della dinamica integrata (lo era per costruzione
Step 4, quando dinamica di derivazione e dinamica integrata
coincidevano). Resta pero' un profilo di guida valido e con senso
fisico (e' esattamente cosi' che viene trattata anche nella pratica
reale: il PEG, derivato dallo stesso principio, viene applicato a
dinamiche via via piu' accurate del modello linearizzato usato per
derivarlo). Il problema inverso (Step 6, ``risolvi_coefficienti_tangente_lineare``)
gia' compensa questo scarto per costruzione, trovando i coefficienti
``A, B`` che centrano il target usando la dinamica REALE (con
``mdot!=0`` e ora con il sollievo centripeto), non quella linearizzata
di derivazione.

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
finale (problema di determinazione di A, B per un target reale: vedi
sezione Step 6 piu' sotto in questo modulo).

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
``tan(theta) = A + B*t`` con ``A``, ``B`` DATI/noti, Step 4) E la
soluzione del problema al contorno per un target di velocita' finale
assegnato (root-finding via shooting, Step 6, sezione piu' sotto). Non
implementa staging durante la fase esoatmosferica (Step 8) ne' la
scelta di ``A, B`` per un target orbitale REALE completo con dati
pubblici (assemblata in ``lanciatore/validazione.py``, Step 7, che
riusa questo modulo cosi' com'e').
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from lanciatore.costanti import R_TERRA
from lanciatore.gravita import accelerazione_gravita


def derivate_stato_tangente_lineare(t, y, A, B, spinta, mdot):
    """Derivate del vettore di stato ``[x, h, vx, vh, m]`` (tangente lineare).

    Equazioni CORRETTE (Ciclo 7, vedi sezione "Correzione fisica" nella
    docstring di modulo per la derivazione completa in coordinate
    polari locali e la dimostrazione di conservazione E/L)::

        dx/dt  = vx
        dh/dt  = vh
        dvx/dt = (spinta/m) * cos(theta(t)) - vh*vx/(R_TERRA+h)
        dvh/dt = (spinta/m) * sin(theta(t)) - g(h) + vx**2/(R_TERRA+h)
        dm/dt  = -mdot

    con ``theta(t) = arctan(A + B*t)`` e ``g(h) =
    lanciatore.gravita.accelerazione_gravita(h)`` (ricalcolata a ogni
    passo, NON piu' congelata a un valore costante -- vedi docstring di
    modulo).

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
        dalle condizioni al contorno sulla velocita' finale (Step 4:
        dati/noti; Step 6: risolti da ``risolvi_coefficienti_tangente_lineare``).
    spinta : float
        Spinta del motore, N.
    mdot : float
        Portata massica, kg/s (positiva; la massa diminuisce). Puo'
        essere 0.0.

    Ritorna
    -------
    numpy.ndarray, shape (5,)
        ``[dx/dt, dh/dt, dvx/dt, dvh/dt, dm/dt]``
    """
    x, h, vx, vh, m = y

    theta = np.arctan(A + B * t)
    a = spinta / m
    r = R_TERRA + h

    dx_dt = vx
    dh_dt = vh
    dvx_dt = a * np.cos(theta) - vh * vx / r
    dvh_dt = a * np.sin(theta) - accelerazione_gravita(h) + vx**2 / r
    dm_dt = -mdot

    return np.array([dx_dt, dh_dt, dvx_dt, dvh_dt, dm_dt])


def integra_tangente_lineare(
    m0, spinta, mdot, A, B, tf, x0=0.0, h0=0.0, vx0=0.0, vh0=0.0
):
    """Integra la dinamica a tangente lineare su un orizzonte temporale fisso.

    Nessun evento terminale: ``tf`` e' fisso e dato (coerente col
    confine di scope di questo modulo -- non e' ancora presente ne'
    staging ne' una condizione di arresto legata al raggiungimento di
    un target orbitale, che appartiene a step successivi).

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
    args = (A, B, spinta, mdot)

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
# a(t) = spinta/m(t) NON e' costante). Con la correzione fisica del Ciclo 7
# non esiste piu' forma chiusa nemmeno nel caso limite mdot=0 (dvx/dt
# dipende ora anch'esso da h/vh, non solo dvh/dt da g) -- vedi docstring di
# modulo, sezione "Correzione fisica".
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
    tf,
    vx_target,
    vh_target,
    x0,
    h0,
    vx0,
    vh0,
):
    """Funzione residuo per fsolve: scarto tra velocita' terminale ottenuta e target.

    Integra la dinamica ESISTENTE (``integra_tangente_lineare``, riusata
    senza alcuna modifica strutturale oltre alla correzione fisica gia'
    applicata a monte) con i coefficienti ``A, B`` correnti e restituisce
    ``[vx(tf) - vx_target, vh(tf) - vh_target]``. Nessuna reimplementazione
    della fisica: questa funzione e' puro glue code tra ``fsolve`` e
    l'integratore.

    Nota (scoperta empiricamente in questo ciclo, non prevista dal piano
    originale): dopo la correzione fisica del Ciclo 7, ``g(h)`` e'
    ricalcolata a ogni passo (non piu' congelata), quindi
    ``accelerazione_gravita`` puo' sollevare ``ValueError`` per ``h<0``
    durante l'integrazione. Con ``g_costante`` (Step 4/6 originali) questo
    non poteva mai accadere (nessun vincolo di dominio su ``h``). Durante
    l'esplorazione dello jacobiano/dei passi di Newton di ``fsolve``,
    specialmente quando il target e' irraggiungibile (nessuna soluzione
    esiste, fsolve tenta passi via via piu' estremi), alcune coppie
    ``(A, B)`` di tentativo producono traiettorie che scendono sotto quota
    zero, facendo sollevare ``ValueError`` -- verificato che questo
    accade per QUALUNQUE target irraggiungibile provato (non e' un
    problema risolvibile scegliendo un margine di quota iniziale diverso,
    e' strutturale: un residuo enorme genera comunque passi di Newton
    enormi). ``fsolve`` non sa gestire eccezioni Python: qui la
    traiettoria fisicamente non valida viene trattata come un punto del
    dominio ``(A, B)`` con residuo enorme ma FINITO, evitando che l'intera
    chiamata crashi con un'eccezione Python grezza.

    Verificato empiricamente (critic-ingegnere, Ciclo 7) che il residuo
    costante NON fornisce alcun gradiente utile a ``fsolve`` dentro la
    regione non valida (lo Jacobiano a differenze finite risulta
    ESATTAMENTE zero li' dentro): il fallback e' un MURO che impedisce il
    crash e trasforma il fallimento in un ``ier`` MINPACK riconoscibile,
    non un PENDIO che guida l'ottimizzatore fuori dalla zona non valida.
    La rete di sicurezza reale contro un esito sbagliato resta il
    controllo combinato gia' presente in ``risolvi_coefficienti_tangente_lineare``
    (flag ``ier`` + residuo effettivo + accordo tra guess primario e
    perturbato): se il target e' irraggiungibile o il guess cade nella
    zona non valida, il risultato finale resta un fallimento esplicito
    (``RuntimeError``), mai una convergenza silenziosa sbagliata
    (verificato su 30 guess di stress test, 0 convergenze errate).
    """
    A, B = coefficienti
    try:
        risultato = integra_tangente_lineare(
            m0, spinta, mdot, A, B, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0
        )
    except ValueError:
        return np.array([1e6, 1e6])
    vx_finale = risultato.y[2, -1]
    vh_finale = risultato.y[3, -1]
    return np.array([vx_finale - vx_target, vh_finale - vh_target])


def risolvi_coefficienti_tangente_lineare(
    m0,
    spinta,
    mdot,
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
    ``integra_tangente_lineare``) su ``[0, tf]`` si ottenga
    ``vx(tf) ~= vx_target`` e ``vh(tf) ~= vh_target``. Con ``mdot != 0``
    reale l'accelerazione ``a(t) = spinta/m(t)`` non e' costante, e con la
    correzione fisica del Ciclo 7 (sollievo centripeto) non esiste forma
    chiusa nemmeno nel caso limite ``mdot = 0`` (vedi docstring di modulo):
    serve necessariamente root-finding numerico in entrambi i casi.

    Fenomeno delle radici multiple/spurie (scoperto empiricamente allo
    Step 6, documentato qui invece di essere nascosto -- vedi STATUS.md,
    Ciclo 6, addendum; ri-verificato numericamente con la fisica corretta
    del Ciclo 7, vedi STATUS.md Ciclo 7 e test dedicato)
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
    ``fsolve``, ma un risultato fisicamente sbagliato. Un tentativo di
    costruire un guess automatico risolvendo la forma chiusa dello Step 4
    (caso ``mdot=0``, ORA NON PIU' DISPONIBILE dopo la correzione del
    Ciclo 7) non risolveva comunque il problema: il problema non e' "quanto
    e' vicino il guess" ma "in quale bacino di attrazione (segno di B) si
    parte" -- informazione che una formula chiusa approssimata non puo'
    fornire in generale senza gia' conoscere la risposta.

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
            args=(m0, spinta, mdot, tf, vx_target, vh_target, x0, h0, vx0, vh0),
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
