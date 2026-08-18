"""Caso di validazione con dati reali (Step 7, Ciclo 7, Fase B).

Assembla la traiettoria completa Stadio 1 (gravity turn, Step 3, riuso non
modificato) -> Stadio 2 (tangente lineare, Step 4/6, equazioni corrette al
Ciclo 7) usando i dati pubblici di un lanciatore reale, e verifica il
delta-v totale contro il benchmark dichiarato in CLAUDE.md (9.1-10.0 km/s
per un'orbita LEO).

Dati veicolo: Falcon 9 Block 5, configurazione ESPENDIBILE
--------------------------------------------------------------------------
Fonte: Wikipedia, "Falcon 9", consultata 2026-08-16,
https://en.wikipedia.org/wiki/Falcon_9 -- tabella delle specifiche
"Falcon 9 Block 5", colonna/riga della variante ESPENDIBILE (non la
variante con recupero/riutilizzo del primo stadio: in configurazione
espendibile il primo stadio non porta il propellente di riserva per
l'atterraggio, quindi tutta la massa di propellente dichiarata e'
effettivamente disponibile per l'ascesa, coerente con l'assenza in questo
progetto di qualunque logica di rientro/atterraggio). Tutti i numeri
sotto (m_prop, m_strut, spinta, Isp di entrambi gli stadi, payload,
diametro) provengono dalla STESSA riga/tabella di quella fonte -- non
mescolati con altre varianti (verificato esplicitamente, non assunto, come
richiesto da STATUS.md Ciclo 7 addendum 3bis punto 2).

Perche' Isp/spinta da VUOTO per ENTRAMBI gli stadi (non SL per lo Stadio 1)
--------------------------------------------------------------------------
Con Isp/spinta a livello del mare (SL) per lo Stadio 1, il delta-v ideale
totale risulta ~8764 m/s, SOTTO il range 9.1-10.0 km/s del benchmark
CLAUDE.md. Con Isp/spinta da vuoto, ~9138 m/s, dentro il range. La scelta
di usare i valori da vuoto NON e' motivata "perche' torna" (STATUS.md
Ciclo 7 addendum 3bis, punto 1, ha richiesto esplicitamente una
giustificazione indipendente dal risultato): con questi dati (spinta
vuoto, kick standard v_kick=50 m/s, kick_angle=2 gradi) lo Stadio 1
brucia fino a una quota di ~79-80 km (verificato numericamente in questo
modulo, vedi ``assembla_traiettoria``), quota alla quale la pressione
atmosferica e' gia' una frazione minima di quella al livello del mare per
la maggior parte della durata del bruciamento -- le prestazioni da vuoto
sono percio' una media pesata-nel-tempo piu' rappresentativa della
prestazione SL per un modello a Isp COSTANTE (semplificazione dichiarata
in CLAUDE.md: nessuna variazione di Isp con la quota e' implementata).
Il numero con Isp1 SL (~8764 m/s, vedi ``ISP_1_SL`` e
``calcola_delta_v_ideale``) resta comunque disponibile nel risultato di
``esegui_validazione`` come limite inferiore di sensibilita', non
nascosto.

Confine delle fasi di guida (STATUS.md Ciclo 7, punto 2 del piano)
--------------------------------------------------------------------------
Stadio 1 = gravity turn INTERO (``guida.integra_gravity_turn``, riuso non
modificato, ``v_kick=50``, ``kick_angle_deg=2.0``); Stadio 2 = tangente
lineare INTERO (equazioni corrette del Ciclo 7 in
``guida_esoatmosferica.py``), ``tf`` = bruciamento completo dello Stadio 2
(``m_prop_2 / mdot_2``). Nessuno staging durante la fase di tangente
lineare (Step 8, esplicitamente fuori scope qui).

Strategia per i coefficienti A, B dello Stadio 2 (STATUS.md Ciclo 7,
addendum 3bis, punto 5 -- DECISIONE DI DESIGN, non un dettaglio)
--------------------------------------------------------------------------
Verificato (dall'orchestratore, prima di questo modulo, e riverificato qui
con l'implementazione reale) che il target esatto (velocita' orbitale a
200 km, componente verticale nulla) NON e' raggiungibile entro la
tolleranza stretta (1e-6) di
``guida_esoatmosferica.risolvi_coefficienti_tangente_lineare`` con
``tf`` fisso: quella funzione (pensata per Step 6, dove la convergenza
esatta era garantita per costruzione del caso di test) qui solleverebbe
sistematicamente ``RuntimeError``. Si usa percio' ``scipy.optimize.minimize``
(Nelder-Mead) sulla somma dei quadrati degli scarti di velocita' finale,
partendo da un guess nella regione indicata dall'addendum
(``A0=-0.04, B0=0.0009``) e verificando la robustezza con un secondo
guess perturbato (stesso principio del doppio guess di Step 6, qui usato
come diagnostica riportata nel risultato, non come gate bloccante: con
Nelder-Mead su un problema di minimi quadrati non c'e' il fenomeno delle
radici multiple/spurie di fsolve, ma la coerenza tra guess diversi resta
un buon controllo di plausibilita'). Lo scarto residuo finale (~30-40 m/s,
~0.4-0.5% del target, verificato empiricamente con guess multipli e
indipendenti che convergono tutti sullo STESSO punto) e' riportato
ESPLICITAMENTE nel risultato, non forzato a zero: e' un limite onesto e
atteso del modello (2 gradi di liberta' A, B per centrare un bersaglio a
2 componenti mentre l'altitudine finale resta libera, con ``tf`` fisso),
non un errore da nascondere.

Scomposizione delle perdite (STATUS.md Ciclo 7, addendum 3bis, punto 3)
--------------------------------------------------------------------------
Identita' (derivata dal reviewer, verificata algebricamente
dall'orchestratore -- i termini centripeti si cancellano esattamente nel
MODULO della velocita', quindi non compaiono qui)::

    dv/dt = a(t)*cos(theta(t)-gamma(t)) - g(h)*sin(gamma(t))

con ``v = sqrt(vx**2+vh**2)``, ``gamma = atan2(vh, vx)``, ``a(t) =
spinta/m(t)`` (accelerazione di SOLA spinta, il drag e' un termine
ADDITIVO separato ``-D(t)/m(t)`` che si somma solo in fase atmosferica,
dove agisce sempre esattamente in direzione ``-v`` per costruzione, quindi
non richiede un fattore coseno separato). Integrando su tutta la
traiettoria::

    v_finale - v_iniziale = INT[a*cos(theta-gamma)]dt - INT[D/m]dt - INT[g*sin(gamma)]dt
                            = (INT[a]dt - INT[a*(1-cos(theta-gamma))]dt) - INT[D/m]dt - INT[g*sin(gamma)]dt
                            = dv_ideale - perdita_manovra - perdita_drag - perdita_gravita

dove ``INT[a]dt = Isp*G0*ln(m0/mf)`` per un bruciamento a portata massica
costante (identita' di Tsiolkovsky, verificata separatamente con
``math.log`` in ``calcola_delta_v_ideale``). Nella fase di gravity turn
(Stadio 1) la spinta e' SEMPRE allineata alla velocita' per ipotesi
(angolo di attacco nullo, Step 3): ``theta = gamma`` identicamente per
costruzione, quindi il contributo alla "perdita di manovra" e'
ESATTAMENTE zero li' (nessun calcolo indipendente necessario: e' zero per
la stessa ragione fisica per cui il gravity turn e' definito cosi'). Nella
fase di tangente lineare (Stadio 2) ``theta = arctan(A+B*t)`` e' noto
analiticamente, ``gamma`` e' ricostruito da ``atan2(vh, vx)`` sullo stato
integrato: qui il contributo NON e' zero (spinta non allineata alla
velocita', guida attiva vera) e viene calcolato con l'integrale DIRETTO
``INT[a*(1-cos(theta-gamma))]dt``, mai per differenza/residuo (vincolo
esplicito del piano).

Nota tecnica importante sulla precisione del check di chiusura
--------------------------------------------------------------------------
Gli oggetti risultato di ``guida.integra_gravity_turn`` e
``guida_esoatmosferica.integra_tangente_lineare`` sono prodotti con
``dense_output=False`` e passo adattivo: i punti restituiti (es. 5 punti
per la Fase A verticale, ~30 per il gravity turn) sono numericamente
ACCURATI in quei punti (coerenti con ``rtol=1e-8`` dell'integratore), ma
sono troppo RADI per una regola del trapezio accurata sugli integrandi
NON lineari di questa scomposizione (``g(h)*sin(gamma)``,
``a*(1-cos(theta-gamma))``): verificato empiricamente che l'identita' di
chiusura, calcolata direttamente sui punti radi, ha uno scarto di ~3.3 m/s
su ~7755 m/s (~0.04%) -- troppo per un "check di precisione
dell'integratore". ``scompone_perdite`` esegue percio' un RE-INTEGRAZIONE
DENSA indipendente (stesse funzioni derivate, stessi parametri/condizioni
iniziali/tf della soluzione gia' calcolata -- nessuna fisica
reimplementata, solo ``solve_ivp`` richiamato di nuovo con
``dense_output=True`` e ``rtol=1e-10``, poi valutato su una griglia fine
di migliaia di punti) usata SOLO per gli integrali di questa
scomposizione, non per il risultato principale della traiettoria (target
raggiunto/scarto, che restano quelli della soluzione originale). Con
questo accorgimento lo scarto di chiusura scende a ~1.5e-4 m/s (~2e-8
relativo), verificato empiricamente -- vedi tolleranza nel test dedicato.
"""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from lanciatore import dinamica, guida
from lanciatore import guida_esoatmosferica as ge
from lanciatore.atmosfera import densita
from lanciatore.costanti import CD, G0, MU_TERRA, R_TERRA
from lanciatore.gravita import accelerazione_gravita

# ---------------------------------------------------------------------------
# Dati veicolo: Falcon 9 Block 5, configurazione ESPENDIBILE.
# Fonte (per OGNI numero sotto): Wikipedia, "Falcon 9", consultata
# 2026-08-16, https://en.wikipedia.org/wiki/Falcon_9 -- tabella delle
# specifiche, variante espendibile (vedi docstring di modulo per la
# motivazione della scelta "espendibile" invece di "riutilizzabile").
# ---------------------------------------------------------------------------

# --- Stadio 1 ---
M_STRUT_1 = 25_600.0            # kg, massa a vuoto (strutturale) del primo stadio
M_PROP_1 = 395_700.0            # kg, massa di propellente del primo stadio
SPINTA_1_VUOTO = 8_227_000.0    # N, spinta totale (9 motori Merlin 1D), VUOTO
ISP_1_VUOTO = 312.0             # s, impulso specifico del primo stadio, VUOTO

# Isp del primo stadio a livello del mare (SL) -- SOLO per il confronto di
# sensibilita' (vedi docstring di modulo), NON usata per la traiettoria
# principale (quella usa sempre ISP_1_VUOTO).
ISP_1_SL = 283.0                # s, impulso specifico del primo stadio, livello del mare

# --- Stadio 2 ---
M_STRUT_2 = 3_900.0             # kg, massa a vuoto (strutturale) del secondo stadio
M_PROP_2 = 92_670.0             # kg, massa di propellente del secondo stadio
SPINTA_2 = 981_000.0            # N, spinta (1 motore Merlin 1D Vacuum)
ISP_2 = 348.0                   # s, impulso specifico del secondo stadio (sempre in vuoto)

# --- Payload e geometria ---
PAYLOAD = 22_800.0              # kg, capacita' massima dichiarata verso LEO, configurazione espendibile
DIAMETRO = 3.7                  # m, diametro del corpo del lanciatore
AREA_RIFERIMENTO = math.pi * (DIAMETRO / 2.0) ** 2  # m^2, area della sezione trasversale (~10.75 m^2)

# ---------------------------------------------------------------------------
# Masse nominali derivate (bookkeeping, STATUS.md Ciclo 7 punto 2):
#   m0 = 540 670 kg, m_vuoto_1 = 144 970 kg, m_ignizione_2 = 119 370 kg,
#   m_vuoto_2 = 26 700 kg (= payload + m_strut_2, coerente per costruzione).
# ---------------------------------------------------------------------------
M0_TOTALE = PAYLOAD + M_STRUT_2 + M_PROP_2 + M_STRUT_1 + M_PROP_1
M_VUOTO_NOMINALE_1 = M0_TOTALE - M_PROP_1  # soglia nominale evento fine_propellente Stadio 1

# ---------------------------------------------------------------------------
# Parametri della legge di guida.
# ---------------------------------------------------------------------------
V_KICK = 50.0            # m/s, riuso Step 3 (kick angle applicato a questa velocita')
KICK_ANGLE_DEG = 2.0     # gradi, riuso Step 3

H_TARGET = 200_000.0     # m, quota orbitale target (LEO circolare, 200 km)

# Guess iniziale per i coefficienti A, B dello Stadio 2 (regione indicata
# dall'addendum 3bis punto 5). Il punto effettivamente trovato dalla
# minimizzazione (~A=-0.00355, B=0.000707) e' MOLTO diverso da questo
# guess, ma verificato che guess diversi (inclusi punti fisicamente
# motivati come A=tan(90-gamma_fine_stadio1)) convergono tutti sullo
# STESSO punto: e' un minimo robusto, non un artefatto del guess.
A0_GUESS_STADIO2 = -0.04
B0_GUESS_STADIO2 = 0.0009

# Guess iniziale per i coefficienti A, B dello Stadio 2 col NUOVO
# obiettivo di quota+verticalita' (STATUS.md Ciclo 10, indagine numerica
# pre-piano): la regione precedente (A0_GUESS_STADIO2 sopra) era tarata
# per il vecchio obiettivo (vx/vh target) e non e' un buon punto di
# partenza per il nuovo. Verificato numericamente che (A0=0.3, B0=0.0)
# converge in modo pulito e a precisione di macchina (fun~=2.4e-14) per
# lo stato di handoff Stadio1->Stadio2 di questo modulo; altri guess
# provati falliscono per ValueError di dominio durante l'esplorazione
# dell'ottimizzatore (comportamento gia' noto, gestito dal fallback
# residuo-costante, vedi _obiettivo_residuo_quota_verticale).
A0_GUESS_STADIO2_QUOTA = 0.3
B0_GUESS_STADIO2_QUOTA = 0.0

# Bonus di velocita' orizzontale dovuto alla rotazione terrestre per un
# sito di lancio equatoriale/subtropicale come Cape Canaveral (28.5°N),
# GIA' quantificato nel disclaimer di CLAUDE.md (sezione aggiunta
# 2026-08-16): 465.1 m/s (velocita' equatoriale di rotazione terrestre)
# proiettata sulla latitudine del sito, cos(28.5°). Il modello di questo
# progetto esclude per costruzione la rotazione terrestre (CLAUDE.md,
# riga "Dinamica"): questa costante e' usata SOLO come termine di
# sensibilita' esplicito in ``verifica_orbita_stabile`` (parametro
# ``bonus_rotazione``), mai applicata di default alla traiettoria
# simulata.
BONUS_ROTAZIONE_TERRESTRE = 465.1 * math.cos(math.radians(28.5))


def converti_v_gamma_a_vx_vh(v, gamma):
    """Converte lo stato polare della Fase B gravity turn in componenti cartesiane.

    Usata alla continuita' Stadio1->Stadio2 (STATUS.md Ciclo 7, addendum
    3bis, punto 4): il gravity turn (Step 3) porta lo stato in forma
    ``[x, h, v, gamma, m]`` (modulo velocita' + flight-path angle), mentre
    la tangente lineare (Step 4/6) richiede ``[x, h, vx, vh, m]``
    (componenti cartesiane). Funzione isolata e testabile separatamente
    dal resto dell'assemblaggio, come richiesto esplicitamente dal piano.

    Parametri
    ---------
    v : float
        Modulo della velocita', m/s.
    gamma : float
        Flight-path angle, RADIANTI (angolo rispetto all'orizzonte
        locale).

    Ritorna
    -------
    (vx, vh) : tuple of float
        Componente orizzontale e verticale della velocita', m/s.
    """
    vx = v * math.cos(gamma)
    vh = v * math.sin(gamma)
    return vx, vh


def massa_ignizione_stadio_2(massa_effettiva_burnout_stadio_1, m_strut1=M_STRUT_1):
    """Massa all'ignizione dello Stadio 2 (continuita' allo staging).

    Massa EFFETTIVA di fine Stadio 1 (dallo stato integrato restituito da
    ``solve_ivp``, NON dalla soglia nominale ``M_VUOTO_NOMINALE_1``) MENO
    la massa strutturale dello Stadio 1 espulsa allo staging -- stesso
    principio gia' stabilito in ``lanciatore/staging.py`` (Step 5) per la
    continuita' tra stadi in fase atmosferica.

    Isolata come funzione a se stante (invece di un'espressione inline
    dentro ``assembla_traiettoria``) apposta per essere testabile in
    isolamento (STATUS.md Ciclo 7, addendum 3bis, punto 4: "test dedicato
    che verifica SEPARATAMENTE... la sottrazione di m_strut_1").

    Parametri
    ---------
    massa_effettiva_burnout_stadio_1 : float
        Massa EFFETTIVA (kg) all'istante dell'evento di fine propellente
        dello Stadio 1, presa da ``risultato.y[4, -1]`` della Fase B del
        gravity turn.
    m_strut1 : float, opzionale
        Massa strutturale dello Stadio 1 espulsa allo staging, kg
        (default ``M_STRUT_1``, dati veicolo di questo modulo).

    Ritorna
    -------
    float
        Massa (kg) all'accensione dello Stadio 2.
    """
    return massa_effettiva_burnout_stadio_1 - m_strut1


def velocita_orbitale_target(h=H_TARGET):
    """Velocita' orbitale circolare alla quota ``h`` (target di iniezione LEO).

    ``v = sqrt(MU_TERRA / (R_TERRA + h))``, stessa formula/costanti gia'
    citate e verificate dallo Step 1 (``lanciatore.costanti``).

    Parametri
    ---------
    h : float, opzionale
        Quota target sul livello del mare, m (default ``H_TARGET`` =
        200 km).

    Ritorna
    -------
    float
        Velocita' orbitale circolare, m/s.
    """
    return math.sqrt(MU_TERRA / (R_TERRA + h))


def calcola_delta_v_ideale(isp1=ISP_1_VUOTO):
    """Delta-v ideale totale (Tsiolkovsky), calcolato con ``math.log`` sulle masse REALI.

    Somma dei due stadi, ciascuno con la propria formula di Tsiolkovsky
    ``Isp*G0*ln(m_inizio/m_fine)`` sulle masse NOMINALI di progetto (non
    sulle masse effettive della traiettoria simulata -- le due dovrebbero
    coincidere a precisione di macchina, dato che la soglia dell'evento
    "fine propellente" e' definita esattamente su queste masse nominali,
    verificato in ``tests/test_validazione.py``).

    Parametri
    ---------
    isp1 : float, opzionale
        Impulso specifico del primo stadio, s (default ``ISP_1_VUOTO``).
        Passare ``ISP_1_SL`` per il confronto di sensibilita' (vedi
        docstring di modulo).

    Ritorna
    -------
    dict
        ``{"dv1": ..., "dv2": ..., "dv_totale": ..., "m0": ..., "m_vuoto_1": ...,
        "m_ignizione_2": ..., "m_vuoto_2": ...}`` (tutte le masse nominali
        intermedie, per trasparenza del bookkeeping).
    """
    m0 = M0_TOTALE
    m_vuoto_1 = M_VUOTO_NOMINALE_1
    m_ignizione_2 = massa_ignizione_stadio_2(m_vuoto_1)  # bookkeeping nominale
    m_vuoto_2 = m_ignizione_2 - M_PROP_2

    dv1 = isp1 * G0 * math.log(m0 / m_vuoto_1)
    dv2 = ISP_2 * G0 * math.log(m_ignizione_2 / m_vuoto_2)

    return {
        "dv1": dv1,
        "dv2": dv2,
        "dv_totale": dv1 + dv2,
        "m0": m0,
        "m_vuoto_1": m_vuoto_1,
        "m_ignizione_2": m_ignizione_2,
        "m_vuoto_2": m_vuoto_2,
    }


def _obiettivo_residuo_quadratico(
    coefficienti, m0, spinta, mdot, tf, x0, h0, vx0, vh0, vx_target, vh_target
):
    """Somma dei quadrati degli scarti di velocita' finale (obiettivo di ``minimize``).

    Stesso principio di fallback della funzione residuo di
    ``guida_esoatmosferica._residuo_velocita_finale`` (Step 6/Ciclo 7):
    ``accelerazione_gravita`` puo' sollevare ``ValueError`` per traiettorie
    di tentativo che scendono sotto quota zero durante l'esplorazione
    dell'ottimizzatore. ``scipy.optimize.minimize`` (a differenza di
    ``fsolve``) non richiede un valore vettoriale ma uno scalare: qui il
    fallback e' un valore scalare enorme ma finito (muro, non pendio,
    stesso comportamento verificato in Ciclo 7 per il caso fsolve).
    """
    A, B = coefficienti
    try:
        risultato = ge.integra_tangente_lineare(
            m0, spinta, mdot, A, B, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0
        )
    except ValueError:
        return 1e12
    if not risultato.success:
        return 1e12
    vx_finale = risultato.y[2, -1]
    vh_finale = risultato.y[3, -1]
    return (vx_finale - vx_target) ** 2 + (vh_finale - vh_target) ** 2


def risolvi_stadio_2_minimi_quadrati(
    m0, spinta, mdot, tf, x0, h0, vx0, vh0, vx_target, vh_target, A0, B0
):
    """Trova A, B dello Stadio 2 per minimi quadrati (non root-finding esatto).

    Vedi docstring di modulo, sezione "Strategia per i coefficienti A, B
    dello Stadio 2", per la motivazione di questa scelta di design
    (target non raggiungibile entro la tolleranza stretta del
    root-finder di Step 6 con questo ``tf`` fisso).

    Esegue la minimizzazione (Nelder-Mead) dal guess primario
    ``(A0, B0)`` E da un secondo guess perturbato, riportando se i due
    concordano (``coerenza_multistart``) -- diagnostica di robustezza,
    NON un gate bloccante: a differenza del root-finding di Step 6
    (``fsolve``, dove un disaccordo tra guess sollevava ``RuntimeError``
    per il fenomeno delle radici spurie), qui un problema di minimi
    quadrati con ``Nelder-Mead`` non ha lo stesso fenomeno di bacini di
    attrazione multipli con radici ugualmente valide: un disaccordo
    verrebbe riportato esplicitamente nel risultato, non nascosto, ma non
    interrompe l'esecuzione (coerente con la decisione di design
    dell'addendum: "si riporta onestamente lo scarto residuo... non si
    forza una convergenza artificiale").

    Ritorna
    -------
    dict
        ``{"A": ..., "B": ..., "risultato": <OdeResult>, "residuo_vx": ...,
        "residuo_vh": ..., "residuo_modulo": ..., "coerenza_multistart": ...,
        "A_secondario": ..., "B_secondario": ...}``.
    """
    args = (m0, spinta, mdot, tf, x0, h0, vx0, vh0, vx_target, vh_target)
    opzioni = {"xatol": 1e-10, "fatol": 1e-12, "maxiter": 20_000, "maxfev": 20_000}

    ris_primario = minimize(
        _obiettivo_residuo_quadratico, x0=[A0, B0], args=args, method="Nelder-Mead", options=opzioni
    )
    A_p, B_p = ris_primario.x

    # Guess perturbato per la diagnostica di robustezza (vedi docstring).
    A0_pert = A0 * 1.3 if A0 != 0.0 else 0.01
    B0_pert = B0 * 1.3 if B0 != 0.0 else 0.0001
    ris_secondario = minimize(
        _obiettivo_residuo_quadratico, x0=[A0_pert, B0_pert], args=args, method="Nelder-Mead", options=opzioni
    )
    A_s, B_s = ris_secondario.x

    coerenza_multistart = bool(np.allclose([A_p, B_p], [A_s, B_s], rtol=1e-2, atol=1e-4))

    risultato = ge.integra_tangente_lineare(m0, spinta, mdot, A_p, B_p, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0)
    vx_finale = risultato.y[2, -1]
    vh_finale = risultato.y[3, -1]

    return {
        "A": float(A_p),
        "B": float(B_p),
        "risultato": risultato,
        "residuo_vx": float(vx_finale - vx_target),
        "residuo_vh": float(vh_finale - vh_target),
        "residuo_modulo": float(math.hypot(vx_finale - vx_target, vh_finale - vh_target)),
        "coerenza_multistart": coerenza_multistart,
        "A_secondario": float(A_s),
        "B_secondario": float(B_s),
    }


def _obiettivo_residuo_quota_verticale(
    coefficienti, m0, spinta, mdot, tf, x0, h0, vx0, vh0, h_target, vh_target
):
    """Obiettivo NORMALIZZATO su quota finale e velocita' verticale finale (Ciclo 10).

    Vedi STATUS.md, Ciclo 10, addendum (reviewer punto 6): a differenza di
    ``_obiettivo_residuo_quadratico`` (che vincola solo vx/vh, lasciando la
    quota finale libera), questo obiettivo vincola ESPLICITAMENTE la quota
    di iniezione (``h_target``) e la componente verticale della velocita'
    (tipicamente ``vh_target=0.0``, iniezione circolare) -- la coppia di
    vincoli che identifica davvero il perigeo/apogeo dell'orbita
    risultante (vedi ``verifica_orbita_stabile``), a differenza del
    vecchio obiettivo che lasciava la quota libera.

    Obiettivo normalizzato (reviewer punto 6: metri^2 e (m/s)^2 sommati
    direttamente sarebbero mal condizionati)::

        ((h(tf) - h_target) / h_target)**2 + (vh(tf) / v_target_circolare)**2

    con ``v_target_circolare = sqrt(MU_TERRA / (R_TERRA + h_target))``
    (``velocita_orbitale_target``, gia' definita in questo modulo) usata
    per dare scala fisicamente sensata al termine di velocita'.

    Stesso principio di fallback di ``_obiettivo_residuo_quadratico``:
    ``ValueError`` di dominio (traiettorie di tentativo che scendono sotto
    quota zero) o integrazione non riuscita -> valore scalare enorme ma
    finito (muro, non pendio).
    """
    A, B = coefficienti
    try:
        risultato = ge.integra_tangente_lineare(
            m0, spinta, mdot, A, B, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0
        )
    except ValueError:
        return 1e12
    if not risultato.success:
        return 1e12
    h_finale = risultato.y[1, -1]
    vh_finale = risultato.y[3, -1]
    v_target_circolare = velocita_orbitale_target(h_target)
    termine_h = ((h_finale - h_target) / h_target) ** 2
    termine_vh = (vh_finale / v_target_circolare) ** 2
    return termine_h + termine_vh


def risolvi_stadio_2_target_quota(
    m0, spinta, mdot, tf, x0, h0, vx0, vh0, h_target, vh_target, A0, B0
):
    """Trova A, B dello Stadio 2 vincolando QUOTA finale e velocita' verticale finale.

    Vedi STATUS.md, Ciclo 10, addendum ("Design finale consolidato"): a
    differenza di ``risolvi_stadio_2_minimi_quadrati`` (target vx/vh,
    quota finale libera -- funzione ESISTENTE, lasciata invariata sopra e
    ancora disponibile per confronto/sensibilita', vedi suo test
    dedicato), questa funzione vincola quota+verticalita', la coppia di
    vincoli che identifica davvero il perigeo/apogeo dell'orbita
    risultante.

    Stesso pattern di doppio guess/diagnostica di robustezza di
    ``risolvi_stadio_2_minimi_quadrati`` (Nelder-Mead dal guess primario
    E da un guess perturbato, ``coerenza_multistart`` riportata come
    diagnostica non bloccante, NON un gate come il ``RuntimeError`` di
    ``fsolve`` in Step 6 -- vedi docstring di quella funzione per la
    motivazione completa, identica qui).

    Ritorna
    -------
    dict
        ``{"A": ..., "B": ..., "risultato": <OdeResult>, "residuo_h": ...,
        "residuo_vh": ..., "residuo_obiettivo": ..., "coerenza_multistart": ...,
        "A_secondario": ..., "B_secondario": ...}``.
    """
    args = (m0, spinta, mdot, tf, x0, h0, vx0, vh0, h_target, vh_target)
    opzioni = {"xatol": 1e-10, "fatol": 1e-14, "maxiter": 20_000, "maxfev": 20_000}

    ris_primario = minimize(
        _obiettivo_residuo_quota_verticale, x0=[A0, B0], args=args, method="Nelder-Mead", options=opzioni
    )
    A_p, B_p = ris_primario.x

    # Guess perturbato per la diagnostica di robustezza (vedi docstring),
    # stesso fallback per B0=0.0 gia' usato in risolvi_stadio_2_minimi_quadrati.
    A0_pert = A0 * 1.3 if A0 != 0.0 else 0.01
    B0_pert = B0 * 1.3 if B0 != 0.0 else 0.0001
    ris_secondario = minimize(
        _obiettivo_residuo_quota_verticale, x0=[A0_pert, B0_pert], args=args, method="Nelder-Mead", options=opzioni
    )
    A_s, B_s = ris_secondario.x

    coerenza_multistart = bool(np.allclose([A_p, B_p], [A_s, B_s], rtol=1e-2, atol=1e-4))

    risultato = ge.integra_tangente_lineare(m0, spinta, mdot, A_p, B_p, tf, x0=x0, h0=h0, vx0=vx0, vh0=vh0)
    h_finale = risultato.y[1, -1]
    vh_finale = risultato.y[3, -1]

    return {
        "A": float(A_p),
        "B": float(B_p),
        "risultato": risultato,
        "residuo_h": float(h_finale - h_target),
        "residuo_vh": float(vh_finale - vh_target),
        "residuo_obiettivo": float(ris_primario.fun),
        "coerenza_multistart": coerenza_multistart,
        "A_secondario": float(A_s),
        "B_secondario": float(B_s),
    }


def verifica_orbita_stabile(h, vx, vh, bonus_rotazione=0.0):
    """Verifica se lo stato finale (h, vx, vh) descrive un'orbita ellittica stabile.

    Meccanica orbitale standard per un campo centrale (energia specifica
    + momento angolare specifico -> semiasse maggiore, eccentricita',
    perigeo, apogeo). Fonte: Curtis, H.D., *Orbital Mechanics for
    Engineering Students* (gia' citato dal progetto dallo Step 1 per
    ``MU_TERRA``/``R_TERRA``), capitoli su energia dell'orbita e momento
    angolare per il problema dei due corpi::

        r = R_TERRA + h
        vx_eff = vx + bonus_rotazione                # bonus SOLO sulla componente orizzontale
        L = r * vx_eff                                # momento angolare specifico
        E = (vx_eff**2 + vh**2) / 2 - MU_TERRA / r    # energia meccanica specifica
        a = -MU_TERRA / (2 * E)                       # semiasse maggiore (orbita ellittica, E<0)
        e = sqrt(1 + 2*E*L**2 / MU_TERRA**2)           # eccentricita'
        r_perigeo = a * (1 - e); r_apogeo = a * (1 + e)

    ``bonus_rotazione`` (default 0.0, il caso rigoroso dichiarato in
    CLAUDE.md: nessuna rotazione terrestre) si applica SOLO a ``vx``
    prima del calcolo -- vedi ``BONUS_ROTAZIONE_TERRESTRE`` in questo
    modulo per il valore usato nel confronto di sensibilita' con
    rotazione terrestre (disclaimer CLAUDE.md, sezione aggiunta
    2026-08-16, rafforzato STATUS.md Ciclo 10).

    Parametri
    ---------
    h : float
        Quota sul livello del mare, m.
    vx, vh : float
        Componenti orizzontale/verticale della velocita', m/s.
    bonus_rotazione : float, opzionale
        Offset costante (m/s) sommato a ``vx`` prima del calcolo (default
        0.0, nessun bonus).

    Ritorna
    -------
    dict
        ``{"semiasse_maggiore": ..., "eccentricita": ..., "perigeo": ...,
        "apogeo": ..., "perigeo_valido": <bool>}``. ``perigeo``/``apogeo``
        sono QUOTE (m) sul livello del mare (``r_perigeo/apogeo - R_TERRA``),
        non distanze geocentriche. ``perigeo_valido = perigeo >= 0``
        (l'orbita risultante non interseca la superficie terrestre).
    """
    r = R_TERRA + h
    vx_eff = vx + bonus_rotazione

    L = r * vx_eff
    E = (vx_eff**2 + vh**2) / 2.0 - MU_TERRA / r
    a = -MU_TERRA / (2.0 * E)
    e = math.sqrt(1.0 + 2.0 * E * L**2 / MU_TERRA**2)

    r_perigeo = a * (1.0 - e)
    r_apogeo = a * (1.0 + e)
    perigeo = r_perigeo - R_TERRA
    apogeo = r_apogeo - R_TERRA

    return {
        "semiasse_maggiore": a,
        "eccentricita": e,
        "perigeo": perigeo,
        "apogeo": apogeo,
        "perigeo_valido": bool(perigeo >= 0.0),
    }


def assembla_traiettoria(A0=A0_GUESS_STADIO2_QUOTA, B0=B0_GUESS_STADIO2_QUOTA, h_target=H_TARGET):
    """Assembla la traiettoria completa Stadio 1 (gravity turn) -> Stadio 2 (tangente lineare).

    Vedi docstring di modulo per il confine delle fasi di guida.
    STATUS.md, Ciclo 10 (design finale consolidato): usa di default
    ``risolvi_stadio_2_target_quota`` (vincola quota finale + velocita'
    verticale finale, invece del vecchio target vx/vh che lasciava la
    quota libera -- vedi Ciclo 9 per la scoperta del problema e Ciclo 10
    per la correzione). ``risolvi_stadio_2_minimi_quadrati`` (Ciclo 7)
    resta nel modulo invariata, non piu' usata qui di default.

    Ritorna
    -------
    dict
        Contiene, tra gli altri campi: ``risultato_1`` (dict grezzo di
        ``guida.integra_gravity_turn``, non post-processato), ``risultato_2``
        (``OdeResult`` grezzo dell'integrazione finale dello Stadio 2 con
        i coefficienti A, B trovati), le masse di continuita' (effettiva
        E nominale, per il confronto esplicito richiesto dal piano),
        ``vx_finale``/``vh_finale``/``v_finale``, il target vx/vh
        INFORMATIVO (``vx_target``/``vh_target``, velocita' circolare
        alla quota target -- NON piu' il vincolo attivo del solutore, solo
        un termine di paragone) con lo scarto residuo corrispondente
        (``residuo_vx``/``residuo_vh``/``residuo_modulo``, informativo,
        atteso NON nullo), il NUOVO residuo di quota attivamente
        vincolato dal solutore (``residuo_h``, atteso vicino a zero), e la
        verifica di stabilita' orbitale dello stato finale sia SENZA
        bonus di rotazione terrestre (``orbita_senza_bonus``, caso
        rigoroso dichiarato in CLAUDE.md) sia CON bonus
        (``orbita_con_bonus``, sensibilita' esplicita, vedi
        ``BONUS_ROTAZIONE_TERRESTRE``).
    """
    mdot1 = SPINTA_1_VUOTO / (ISP_1_VUOTO * G0)

    # --- Stadio 1: gravity turn intero, riuso non modificato di Step 3. ---
    risultato_1 = guida.integra_gravity_turn(
        M0_TOTALE,
        M_VUOTO_NOMINALE_1,
        SPINTA_1_VUOTO,
        ISP_1_VUOTO,
        CD,
        AREA_RIFERIMENTO,
        V_KICK,
        KICK_ANGLE_DEG,
    )
    gt1 = risultato_1["gravity_turn"]

    # Diagnostica non bloccante (stesso principio di staging.py, Step 5:
    # se lo stadio si e' fermato per un evento difensivo invece che per
    # fine-propellente, e' un segnale da riportare al chiamante/test, non
    # da nascondere ne' da far crashare qui dentro).
    stadio_1_fine_propellente = len(gt1.t_events[0]) > 0

    x1, h1, v1, gamma1, m1_effettiva = gt1.y[:, -1]
    vx1, vh1 = converti_v_gamma_a_vx_vh(v1, gamma1)

    # --- Continuita' Stadio1 -> Stadio2 (addendum 3bis, punto 4). ---
    m_ignizione_2_effettiva = massa_ignizione_stadio_2(m1_effettiva)
    m_ignizione_2_nominale = massa_ignizione_stadio_2(M_VUOTO_NOMINALE_1)

    # --- Stadio 2: tangente lineare intero. ---
    mdot2 = SPINTA_2 / (ISP_2 * G0)
    tf2 = M_PROP_2 / mdot2

    # Target vx/vh INFORMATIVO (velocita' circolare alla quota target):
    # non e' piu' il vincolo attivo del solutore (Ciclo 10), ma resta
    # calcolato e riportato per confronto -- vedi ``residuo_vx``/
    # ``residuo_vh``/``residuo_modulo`` sotto, tutti informativi.
    vx_target = velocita_orbitale_target(h_target)
    vh_target = 0.0

    soluzione_2 = risolvi_stadio_2_target_quota(
        m_ignizione_2_effettiva, SPINTA_2, mdot2, tf2, x1, h1, vx1, vh1, h_target, vh_target, A0, B0
    )
    risultato_2 = soluzione_2["risultato"]
    h_finale = risultato_2.y[1, -1]
    vx_finale = risultato_2.y[2, -1]
    vh_finale = risultato_2.y[3, -1]

    # Verifica di stabilita' orbitale dello stato finale (STATUS.md,
    # Ciclo 10): sia SENZA bonus di rotazione terrestre (caso rigoroso
    # dichiarato in CLAUDE.md, ci si aspetta perigeo_valido=False con
    # questi dati) sia CON bonus (sensibilita' esplicita, ci si aspetta
    # perigeo_valido=True) -- entrambi riportati, nessuno nascosto.
    orbita_senza_bonus = verifica_orbita_stabile(h_finale, vx_finale, vh_finale)
    orbita_con_bonus = verifica_orbita_stabile(
        h_finale, vx_finale, vh_finale, bonus_rotazione=BONUS_ROTAZIONE_TERRESTRE
    )

    return {
        "risultato_1": risultato_1,
        "risultato_2": risultato_2,
        "stadio_1_fine_propellente": stadio_1_fine_propellente,
        "x1": x1,
        "h1": h1,
        "v1": v1,
        "gamma1": gamma1,
        "vx1": vx1,
        "vh1": vh1,
        "m1_effettiva": m1_effettiva,
        "m_ignizione_2_effettiva": m_ignizione_2_effettiva,
        "m_ignizione_2_nominale": m_ignizione_2_nominale,
        "mdot2": mdot2,
        "tf2": tf2,
        "h_target": h_target,
        "A": soluzione_2["A"],
        "B": soluzione_2["B"],
        "coerenza_multistart": soluzione_2["coerenza_multistart"],
        "vx_target": vx_target,
        "vh_target": vh_target,
        "h_finale": float(h_finale),
        "vx_finale": float(vx_finale),
        "vh_finale": float(vh_finale),
        "v_finale": float(math.hypot(vx_finale, vh_finale)),
        "residuo_vx": float(vx_finale - vx_target),
        "residuo_vh": float(vh_finale - vh_target),
        "residuo_modulo": float(math.hypot(vx_finale - vx_target, vh_finale - vh_target)),
        "residuo_h": soluzione_2["residuo_h"],
        "orbita_senza_bonus": orbita_senza_bonus,
        "orbita_con_bonus": orbita_con_bonus,
    }


def scompone_perdite(assemblaggio, n_punti_verticale=5000, n_punti_gravity_turn=20_000, n_punti_tangente=20_000):
    """Scompone il delta-v ideale nei 4 termini (ideale, gravita', drag, manovra).

    Vedi docstring di modulo per l'identita' matematica e per la nota
    tecnica sulla ri-integrazione densa (necessaria per la precisione del
    check di chiusura). Nessuna fisica reimplementata: le tre
    ri-integrazioni usano le STESSE funzioni derivate gia' testate
    (``lanciatore.dinamica.derivate_stato``,
    ``lanciatore.guida.derivate_stato_gravity_turn``,
    ``lanciatore.guida_esoatmosferica.derivate_stato_tangente_lineare``),
    con gli STESSI parametri/condizioni iniziali/tf gia' usati per
    produrre ``assemblaggio`` (nessun nuovo dato fisico introdotto).

    Parametri
    ---------
    assemblaggio : dict
        Risultato di ``assembla_traiettoria``.
    n_punti_verticale, n_punti_gravity_turn, n_punti_tangente : int, opzionali
        Numero di punti della griglia fine di valutazione per ciascun
        segmento (default verificati empiricamente sufficienti a
        raggiungere un errore di chiusura di ~1e-4 m/s, vedi docstring di
        modulo).

    Ritorna
    -------
    dict
        ``{"dv_ideale_totale": ..., "perdita_gravita": ..., "perdita_drag": ...,
        "perdita_manovra": ..., "v_iniziale": ..., "v_finale": ...,
        "residuo_chiusura": ..., "residuo_chiusura_relativo": ...}``.
        ``residuo_chiusura = (dv_ideale_totale - perdita_gravita - perdita_drag
        - perdita_manovra) - (v_finale - v_iniziale)``: deve essere vicino a
        zero entro la precisione dell'integratore (verificato nel test
        dedicato, tolleranza esplicita e stretta).
    """
    risultato_1 = assemblaggio["risultato_1"]
    fa = risultato_1["fase_verticale"]
    gt1 = risultato_1["gravity_turn"]
    risultato_2 = assemblaggio["risultato_2"]

    mdot1 = SPINTA_1_VUOTO / (ISP_1_VUOTO * G0)
    mdot2 = assemblaggio["mdot2"]
    A_opt = assemblaggio["A"]
    B_opt = assemblaggio["B"]

    # --- Segmento 1: ascesa verticale pura (Fase A, Stadio 1). ---
    sol1 = solve_ivp(
        dinamica.derivate_stato,
        t_span=(0.0, fa.t[-1]),
        y0=fa.y[:, 0],
        method="RK45",
        args=(SPINTA_1_VUOTO, ISP_1_VUOTO, CD, AREA_RIFERIMENTO, mdot1, M_VUOTO_NOMINALE_1),
        rtol=1e-10,
        atol=[1e-6, 1e-9, 1e-6],
        dense_output=True,
    )
    t1 = np.linspace(0.0, fa.t[-1], n_punti_verticale)
    h_1, v_1, m_1 = sol1.sol(t1)
    D_1 = 0.5 * densita(h_1) * CD * AREA_RIFERIMENTO * v_1**2
    # Ascesa verticale pura: gamma = 90 gradi per costruzione (nessuna
    # componente orizzontale), spinta verticale allineata alla velocita'
    # verticale => theta = gamma per costruzione => perdita di manovra
    # ESATTAMENTE zero qui (nessun calcolo indipendente necessario).
    grav_1 = accelerazione_gravita(h_1)  # * sin(90 gradi) = 1
    drag_1 = D_1 / m_1

    # --- Segmento 2: gravity turn (Fase B, Stadio 1). ---
    sol2 = solve_ivp(
        guida.derivate_stato_gravity_turn,
        t_span=(0.0, gt1.t[-1]),
        y0=gt1.y[:, 0],
        method="RK45",
        args=(SPINTA_1_VUOTO, ISP_1_VUOTO, CD, AREA_RIFERIMENTO, mdot1),
        rtol=1e-10,
        atol=[1e-6, 1e-6, 1e-9, 1e-12, 1e-6],
        dense_output=True,
    )
    t2 = np.linspace(0.0, gt1.t[-1], n_punti_gravity_turn)
    x_2, h_2, v_2, gamma_2, m_2 = sol2.sol(t2)
    D_2 = 0.5 * densita(h_2) * CD * AREA_RIFERIMENTO * v_2**2
    grav_2 = accelerazione_gravita(h_2) * np.sin(gamma_2)
    drag_2 = D_2 / m_2
    # Gravity turn: spinta sempre allineata alla velocita' per ipotesi
    # (Step 3) => theta = gamma per costruzione => perdita di manovra zero.

    # --- Segmento 3: tangente lineare (Stadio 2). ---
    sol3 = solve_ivp(
        ge.derivate_stato_tangente_lineare,
        t_span=(0.0, assemblaggio["tf2"]),
        y0=risultato_2.y[:, 0],
        method="RK45",
        args=(A_opt, B_opt, SPINTA_2, mdot2),
        rtol=1e-10,
        atol=[1e-6, 1e-6, 1e-9, 1e-9, 1e-6],
        dense_output=True,
    )
    t3 = np.linspace(0.0, assemblaggio["tf2"], n_punti_tangente)
    x_3, h_3, vx_3, vh_3, m_3 = sol3.sol(t3)
    gamma_3 = np.arctan2(vh_3, vx_3)
    theta_3 = np.arctan(A_opt + B_opt * t3)
    a_3 = SPINTA_2 / m_3
    grav_3 = accelerazione_gravita(h_3) * np.sin(gamma_3)
    man_3 = a_3 * (1.0 - np.cos(theta_3 - gamma_3))
    # Nessun drag in fase esoatmosferica (ipotesi di progetto, CLAUDE.md).

    perdita_gravita = float(np.trapezoid(grav_1, t1) + np.trapezoid(grav_2, t2) + np.trapezoid(grav_3, t3))
    perdita_drag = float(np.trapezoid(drag_1, t1) + np.trapezoid(drag_2, t2))
    perdita_manovra = float(np.trapezoid(man_3, t3))  # zero per costruzione nei segmenti 1-2

    dv_ideale_totale = calcola_delta_v_ideale()["dv_totale"]

    v_iniziale = 0.0
    v_finale = assemblaggio["v_finale"]

    lato_sinistro = dv_ideale_totale - perdita_gravita - perdita_drag - perdita_manovra
    lato_destro = v_finale - v_iniziale
    residuo_chiusura = lato_sinistro - lato_destro

    return {
        "dv_ideale_totale": dv_ideale_totale,
        "perdita_gravita": perdita_gravita,
        "perdita_drag": perdita_drag,
        "perdita_manovra": perdita_manovra,
        "v_iniziale": v_iniziale,
        "v_finale": v_finale,
        "residuo_chiusura": residuo_chiusura,
        "residuo_chiusura_relativo": residuo_chiusura / lato_destro,
    }


def esegui_validazione(A0=A0_GUESS_STADIO2_QUOTA, B0=B0_GUESS_STADIO2_QUOTA):
    """Esegue l'intera pipeline di validazione e riporta tutti i numeri rilevanti.

    Nota (pulizia Ciclo 10): il default deve usare i guess tarati per
    l'obiettivo quota/verticalita' (``*_QUOTA``), coerente con
    ``assembla_traiettoria`` che dal Ciclo 10 usa
    ``risolvi_stadio_2_target_quota`` di default. Prima di questa
    correzione i default erano i vecchi guess ``A0_GUESS_STADIO2``/
    ``B0_GUESS_STADIO2`` (tarati per l'obiettivo vx/vh, non piu' quello
    di default) — l'ottimizzatore convergeva comunque allo stesso punto
    (verificato dal critic-ingegnere), quindi non era un bug numerico,
    ma il nome/intento del parametro non corrispondeva al suo uso reale.

    Ritorna
    -------
    dict
        ``{"assemblaggio": ..., "perdite": ..., "delta_v_ideale": ...,
        "delta_v_ideale_sensibilita_isp1_sl": ...}``.
    """
    assemblaggio = assembla_traiettoria(A0=A0, B0=B0)
    perdite = scompone_perdite(assemblaggio)
    delta_v_ideale = calcola_delta_v_ideale()
    delta_v_ideale_sensibilita_isp1_sl = calcola_delta_v_ideale(isp1=ISP_1_SL)

    return {
        "assemblaggio": assemblaggio,
        "perdite": perdite,
        "delta_v_ideale": delta_v_ideale,
        "delta_v_ideale_sensibilita_isp1_sl": delta_v_ideale_sensibilita_isp1_sl,
    }
