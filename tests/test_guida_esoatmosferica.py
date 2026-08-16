"""Verifica della guida a tangente lineare, fase esoatmosferica (Step 4/6,
fisica corretta al Ciclo 7).

Parametri di test marcati "# PROVVISORIO" (vedi STATUS.md, Ciclo 4,
punto 4): non sono dati di un lanciatore reale (riservato allo Step 6/7),
ne' una soluzione del problema al contorno A/B-per-target orbitale reale
(assemblata in lanciatore/validazione.py, fuori scope qui).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo
modulo (riservato a lanciatore/validazione.py).

IMPORTANTE (Ciclo 7): dopo la correzione fisica (sollievo centripeto,
vedi docstring di lanciatore/guida_esoatmosferica.py), NON esiste piu'
alcuna forma chiusa per vx(tf)/vh(tf), nemmeno nel caso limite mdot=0
(dvx/dt ora dipende anch'esso da h/vh, non solo dvh/dt da g). I vecchi
test basati su numpy.arcsinh/numpy.sqrt (Step 4) sono stati SUPERATI, non
adattati con tolleranze piu' larghe: sono stati riscritti secondo la
nuova strategia di verifica descritta in STATUS.md, Ciclo 7, addendum
punto 1bis (check di conservazione nel limite spinta=0, e target sempre
generati per integrazione avanti, mai forma chiusa).
"""

import math

import numpy as np
import pytest

from lanciatore.costanti import G0, MU_TERRA, R_TERRA
from lanciatore.guida_esoatmosferica import (
    derivate_stato_tangente_lineare,
    integra_tangente_lineare,
    risolvi_coefficienti_tangente_lineare,
)

# PROVVISORIO -- Step 4, parametri di test generici (vedi STATUS.md,
# Ciclo 4, punto 4). mdot=0 resta un caso limite utile per isolare la
# dinamica di puro impulso da quella di deplezione di massa, ma NON
# produce piu' una forma chiusa dopo la correzione del Ciclo 7 (vedi
# docstring di modulo sopra): i confronti sotto sono tutti per
# integrazione, non per formula analitica.
A = 0.3
B = -0.002  # 1/s
TF = 20.0  # s
M0 = 1000.0  # kg
SPINTA = 40_000.0  # N -> a0 = spinta/m0 = 40 m/s^2
X0 = H0 = VX0 = VH0 = 0.0


def test_1_verifica_indipendente_integratore_alternativo():
    # Sostituisce il vecchio criterio 1 ("confronto forma chiusa",
    # asinh/sqrt), non piu' disponibile dopo la correzione fisica del
    # Ciclo 7 (vedi docstring di modulo). La nuova verifica e'
    # indipendente da integra_tangente_lineare/solve_ivp nel senso che
    # usa un integratore diverso, scritto qui a mano (RK4 esplicito a
    # passo fisso, molto fine): cattura bug di "wiring" nella funzione
    # di libreria (ordine di args, condizioni iniziali, tolleranze
    # dell'integratore) che un confronto della funzione contro se stessa
    # non potrebbe mai rivelare. Non e' pero' indipendente da un
    # eventuale errore nella FORMULA delle derivate stessa (quello e'
    # coperto, in modo davvero indipendente e non tautologico, dal
    # criterio 2 sotto, che ricalcola g_eff da zero).
    passo = 1e-3  # s, fine abbastanza da garantire convergenza RK4
    n_passi = int(round(TF / passo))

    y = np.array([X0, H0, VX0, VH0, M0])
    t = 0.0
    for _ in range(n_passi):
        k1 = derivate_stato_tangente_lineare(t, y, A, B, SPINTA, 0.0)
        k2 = derivate_stato_tangente_lineare(t + passo / 2, y + passo / 2 * k1, A, B, SPINTA, 0.0)
        k3 = derivate_stato_tangente_lineare(t + passo / 2, y + passo / 2 * k2, A, B, SPINTA, 0.0)
        k4 = derivate_stato_tangente_lineare(t + passo, y + passo * k3, A, B, SPINTA, 0.0)
        y = y + (passo / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += passo

    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    assert risultato.y[2, -1] == pytest.approx(y[2], rel=1e-5)
    assert risultato.y[3, -1] == pytest.approx(y[3], rel=1e-5)


def test_2_identita_algebrica_diretta():
    # Criterio 2: valutando derivate_stato_tangente_lineare a un
    # istante/stato arbitrario, dvx/dt e dvh/dt restituiti devono
    # coincidere con la formula CORRETTA (Ciclo 7, entrambi i termini di
    # sollievo centripeto), ricalcolata QUI in modo indipendente
    # dall'implementazione (g_eff da MU_TERRA/(R_TERRA+h)**2, non
    # importando accelerazione_gravita) -- altrimenti il test sarebbe
    # tautologico (copiare la stessa formula da entrambe le parti).
    t_test = 7.3
    y = np.array([123.0, 456.0, 78.0, 9.0, 850.0])
    spinta_test = 55_000.0
    x, h, vx, vh, m = y

    theta_atteso = np.arctan(A + B * t_test)
    r = R_TERRA + h
    g_atteso = MU_TERRA / r**2  # ricalcolo indipendente, stessa fonte fisica ma non la stessa riga di codice
    dvx_atteso = (spinta_test / m) * np.cos(theta_atteso) - vh * vx / r
    dvh_atteso = (spinta_test / m) * np.sin(theta_atteso) - g_atteso + vx**2 / r

    derivate = derivate_stato_tangente_lineare(
        t_test, y, A=A, B=B, spinta=spinta_test, mdot=12.0
    )

    assert derivate[2] == pytest.approx(dvx_atteso, rel=1e-12)
    assert derivate[3] == pytest.approx(dvh_atteso, rel=1e-12)


def test_3_monotonia_quota_caso_nominale():
    # Criterio 3: h(t) deve essere strettamente crescente nel caso di
    # test nominale. Con la correzione del Ciclo 7 il sollievo
    # centripeto in dvh/dt (+vx**2/(R_TERRA+h)) e' un contributo
    # ulteriormente POSITIVO rispetto alla versione precedente, quindi
    # se il criterio valeva prima (dvh/dt>0 per l'intera durata) resta
    # a fortiori valido ora. Se questo test fallisse sarebbe un errore
    # da investigare, non da nascondere (regola generale CLAUDE.md).
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    h = risultato.y[1]
    assert np.all(np.diff(h) > 0)


def test_4_nessun_nan_inf():
    # Criterio 4: nessun NaN/inf su tutta la traiettoria (caso
    # nominale, mdot=0).
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    for componente in risultato.y:
        assert np.all(np.isfinite(componente))


def test_5_massa_costante_quando_mdot_zero():
    # Criterio 5: verifica banale ma diretta che con mdot=0 la massa
    # resti esattamente m0 per tutta la traiettoria (nessun bug che la
    # fa decrescere comunque).
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    m = risultato.y[4]
    assert m == pytest.approx(M0, rel=0.0, abs=1e-6)


def test_6_caso_con_deplezione_massa_mdot_diverso_zero():
    # Criterio 6 (sanity, non hand-check): con parametri di deplezione
    # realistici (isp=300s coerente con G0, come negli step
    # precedenti), nessun NaN/inf, massa decrescente correttamente. Il
    # confronto analitico esatto non e' piu' disponibile in nessun caso
    # (vedi docstring di modulo) -- questo test resta un controllo di
    # sanita' fisica, non un hand-check numerico.
    isp = 300.0  # s
    mdot = SPINTA / (isp * G0)

    # Protezione obbligatoria (gia' presente dallo Step 4/6): verificare
    # ESPLICITAMENTE a monte che mdot*tf resti ben al di sotto di m0,
    # cosi' che spinta/m non rischi mai una singolarita' durante
    # l'integrazione a tf fisso.
    assert mdot * TF < 0.5 * M0, (
        "mdot*tf non e' sufficientemente piccolo rispetto a m0: rischio di "
        "massa residua troppo bassa entro tf con questi parametri di test."
    )

    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=mdot, A=A, B=B, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    for componente in risultato.y:
        assert np.all(np.isfinite(componente))

    m = risultato.y[4]
    assert np.all(np.diff(m) < 0)  # massa strettamente decrescente
    assert m[-1] == pytest.approx(M0 - mdot * TF, rel=1e-6)


# ---------------------------------------------------------------------------
# Test di conservazione nel limite spinta=0 (Ciclo 7, addendum 1bis, punto
# 7): drag gia' assente in fase esoatmosferica per costruzione, quindi con
# spinta=0 la dinamica e' quella puramente balistica/centrale dimostrata
# nella docstring di modulo (conservazione ESATTA di L ed E). Questi due
# test sostituiscono concettualmente il vecchio confronto in forma chiusa
# come check "primario" di correttezza fisica delle equazioni integrate.
#
# Condizioni iniziali PROVVISORIE ma fisicamente sensate: quota ~150 km e
# velocita' orizzontale ~7000 m/s (vicina al regime orbitale) sono scelte
# apposta per rendere il contributo del sollievo centripeto NON
# trascurabile nel check (a basse velocita' il test sarebbe debole: la
# correzione varrebbe comunque, ma sarebbe una frazione minuscola della
# dinamica totale, quindi un bug nel segno/fattore del nuovo termine
# potrebbe passare inosservato entro le tolleranze dell'integratore).
# Verificato numericamente (script ad-hoc) che la quota resta positiva per
# l'intera durata (non serve evento di sicurezza qui, nessun evento e'
# comunque previsto da questo modulo per costruzione).
# ---------------------------------------------------------------------------

H0_CONS = 150_000.0  # m
VX0_CONS = 7000.0  # m/s
VH0_CONS = 100.0  # m/s
TF_CONS = 200.0  # s


def test_11_conservazione_momento_angolare_spinta_zero():
    # L = vx*(R_TERRA+h) e' un'identita' ESATTA nel limite balistico
    # (dimostrato algebricamente nella docstring di modulo a partire
    # dalle equazioni corrette), quindi la tolleranza accettabile e'
    # quella dell'integratore (rtol=1e-8 in integra_tangente_lineare),
    # non un margine fisico allargato.
    risultato = integra_tangente_lineare(
        M0, spinta=0.0, mdot=0.0, A=0.0, B=0.0, tf=TF_CONS,
        x0=0.0, h0=H0_CONS, vx0=VX0_CONS, vh0=VH0_CONS,
    )
    assert risultato.success

    vx = risultato.y[2]
    h = risultato.y[1]
    L = vx * (R_TERRA + h)

    variazione_relativa = np.max(np.abs(L - L[0]) / np.abs(L[0]))
    assert variazione_relativa < 1e-6


def test_12_conservazione_energia_spinta_zero():
    # E = (vx**2+vh**2)/2 - MU_TERRA/(R_TERRA+h), stessa logica del
    # criterio precedente (identita' esatta, tolleranza dell'integratore
    # non ampliata).
    risultato = integra_tangente_lineare(
        M0, spinta=0.0, mdot=0.0, A=0.0, B=0.0, tf=TF_CONS,
        x0=0.0, h0=H0_CONS, vx0=VX0_CONS, vh0=VH0_CONS,
    )
    assert risultato.success

    vx = risultato.y[2]
    vh = risultato.y[3]
    h = risultato.y[1]
    E = 0.5 * (vx**2 + vh**2) - MU_TERRA / (R_TERRA + h)

    variazione_relativa = np.max(np.abs(E - E[0]) / np.abs(E[0]))
    assert variazione_relativa < 1e-6


# ---------------------------------------------------------------------------
# Step 6 -- problema inverso: risolvi_coefficienti_tangente_lineare.
#
# Parametri PROVVISORI condivisi da questi test (vedi STATUS.md, Ciclo 6 +
# addendum): mdot REALE (isp=300s), A_vero/B_vero riusano deliberatamente
# le stesse costanti A, B definite sopra per continuita' col Ciclo 4.
#
# NOTA (scoperta durante l'implementazione di questo ciclo, non prevista
# dal piano -- vedi anche il report finale): con H0=0.0 (usato dai test
# 1-6 sopra) il problema inverso diventa numericamente fragile con la
# fisica corretta. accelerazione_gravita(h) solleva ValueError per h<0
# (comportamento gia' esistente, invariato), ma ora g(h) e' ricalcolata a
# ogni passo invece di essere congelata: durante l'esplorazione dello
# jacobiano di fsolve, alcune coppie (A, B) vicine al guess (incluso il
# guess raccomandato A0=A_vero*0.7 dall'addendum Ciclo 6) producono
# dvh/dt(0) < 0 con vh0=0/h0=0 (spinta verticale iniziale insufficiente a
# vincere g(0) per QUEL guess specifico, non per la soluzione vera) e la
# traiettoria scende sotto quota zero nei primi istanti, facendo
# sollevare ValueError invece di un semplice residuo "cattivo" per
# fsolve. Con g_costante (Step 4/6 originali) questo non accadeva mai:
# g era un numero fisso, senza vincolo di dominio su h. Per i test del
# problema inverso si usa quindi H0_INV=1000.0 m (un piccolo margine di
# quota positiva, comunque coerente con "fase esoatmosferica": questo
# modulo non ha mai preteso che H0=0 fosse una quota esoatmosferica
# realistica, era solo un placeholder numerico "PROVVISORIO"), verificato
# empiricamente sufficiente a evitare il problema per tutti i guess usati
# in questo file.
# ---------------------------------------------------------------------------

ISP_STEP6 = 300.0  # s
MDOT_STEP6 = SPINTA / (ISP_STEP6 * G0)  # ~= 13.609 kg/s

X0_INV = 0.0
H0_INV = 1000.0  # m -- vedi nota sopra
VX0_INV = 0.0
VH0_INV = 0.0

# Guess raccomandato dall'addendum Ciclo 6: A0 = A_vero*0.7, B0 = B_vero*0.7
# (stesso segno del vero, entro lo stesso bacino di attrazione). Ancora
# valido con la fisica corretta (verificato numericamente in questo ciclo).
A0_GUESS_CORRETTO = A * 0.7
B0_GUESS_CORRETTO = B * 0.7


def test_7_recupero_non_tautologico_mdot_reale():
    # Criterio non tautologico (STATUS.md Ciclo 6, punto 2.4): il target
    # e' generato integrando AVANTI con A_vero, B_vero e mdot REALE (via
    # integra_tangente_lineare) -- risultato NUMERICO, non forma chiusa
    # (non esiste in nessun caso dopo il Ciclo 7). Il solver deve
    # RECUPERARE A_vero, B_vero partendo da un guess diverso ma nello
    # stesso bacino (A0, B0 scalati 0.7x), non solo azzerare il residuo
    # sul target.
    assert MDOT_STEP6 * TF < 0.5 * M0  # margine di massa ampio

    target = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A, B=B, tf=TF,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )
    assert target.success
    vx_target = target.y[2, -1]
    vh_target = target.y[3, -1]

    A_risolto, B_risolto = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, MDOT_STEP6, TF, vx_target, vh_target,
        A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )

    # Recupero dei parametri noti (non solo del target), tolleranza stretta.
    assert A_risolto == pytest.approx(A, rel=1e-4)
    assert B_risolto == pytest.approx(B, rel=1e-4)

    # Verifica aggiuntiva: re-integrando con A, B risolti si riottiene il
    # target entro la tolleranza di residuo del solver.
    riverifica = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A_risolto, B=B_risolto, tf=TF,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )
    assert riverifica.y[2, -1] == pytest.approx(vx_target, abs=1e-4)
    assert riverifica.y[3, -1] == pytest.approx(vh_target, abs=1e-4)


def test_8_companion_mdot_zero_target_da_integrazione_avanti():
    # Test companion (STATUS.md Ciclo 6, addendum, "Test companion non
    # circolare") RIPROGETTATO al Ciclo 7 (addendum 1bis, punto 5): il
    # target NON puo' piu' essere generato con la formula chiusa
    # asinh/sqrt (non esiste piu', vedi docstring di modulo), quindi si
    # genera anch'esso per integrazione avanti, come nel criterio 7, ma
    # nel caso limite mdot=0 (invece di mdot reale) -- unifica la
    # strategia gia' usata per mdot!=0, come indicato esplicitamente
    # nell'addendum.
    #
    # Nota di trasparenza: questo test NON e' piu' "non circolare" nel
    # senso stretto originale (generazione del target indipendente da
    # integra_tangente_lineare) -- quella indipendenza richiedeva la
    # forma chiusa, che la correzione fisica del Ciclo 7 ha eliminato
    # strutturalmente in ogni caso (dvx/dt dipende ora da h/vh anche con
    # mdot=0). Il valore residuo di questo test e' verificare che il
    # solver funzioni correttamente anche nel regime limite mdot=0 (un
    # percorso di codice fisicamente distinto da mdot!=0, es. massa mai
    # decrescente), non piu' l'indipendenza dei dati.
    target = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, tf=TF,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )
    assert target.success
    vx_target = target.y[2, -1]
    vh_target = target.y[3, -1]

    A_risolto, B_risolto = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, 0.0, TF, vx_target, vh_target,
        A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )

    assert A_risolto == pytest.approx(A, rel=1e-4)
    assert B_risolto == pytest.approx(B, rel=1e-4)


def test_9_fallimento_esplicito_target_irraggiungibile():
    # Criterio di fallimento (STATUS.md Ciclo 6, addendum): vh_target ben
    # oltre il limite di Tsiolkovsky per questi parametri deve sollevare
    # RuntimeError, non un risultato silenzioso e sbagliato. Confermato
    # dal reviewer (Ciclo 7, addendum 1bis, punto 6): resta valido solo
    # di firma, nessuna modifica concettuale -- il limite di Tsiolkovsky
    # e' indipendente dal modello di gravita' per costruzione (dipende
    # solo da Isp e rapporto di massa).
    m_finale = M0 - MDOT_STEP6 * TF
    dv_ideale_tsiolkovsky = ISP_STEP6 * G0 * math.log(M0 / m_finale)

    vh_target_irraggiungibile = 5000.0  # m/s
    assert vh_target_irraggiungibile > dv_ideale_tsiolkovsky  # verifica esplicita del limite

    vx_target_qualsiasi = 899.267274  # non rilevante: vh da solo gia' impossibile

    with pytest.raises(RuntimeError):
        risolvi_coefficienti_tangente_lineare(
            M0, SPINTA, MDOT_STEP6, TF,
            vx_target_qualsiasi, vh_target_irraggiungibile,
            A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
            x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
        )


def test_10_documentazione_fenomeno_radici_spurie():
    # Documenta il fenomeno scoperto allo Step 6 (STATUS.md Ciclo 6,
    # addendum): un guess con B0 di segno OPPOSTO al vero converge
    # (residuo piccolo, nessuna RuntimeError) ma su A, B DIVERSI da
    # A_vero, B_vero. RI-VERIFICATO ESPLICITAMENTE in questo ciclo
    # (Ciclo 7, addendum 1bis, punto 8) con la fisica corretta: il
    # fenomeno SI RIPRESENTA (non era garantito a priori che la
    # topologia del sistema residuo restasse la stessa dopo aver
    # aggiunto i due nuovi termini). Non e' un test che deve "fallire"
    # nel senso pytest -- verifica che il fenomeno e' riproducibile, per
    # proteggere contro una futura modifica che lo nasconda senza che
    # nessuno se ne accorga.
    target = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A, B=B, tf=TF,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )
    assert target.success
    vx_target = target.y[2, -1]
    vh_target = target.y[3, -1]

    # Guess con B0 di segno OPPOSTO al vero (B_vero = -0.002).
    A0_guess_sbagliato = A * 0.7
    B0_guess_sbagliato = abs(B) * 0.7  # segno invertito rispetto a B_vero

    A_spurio, B_spurio = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, MDOT_STEP6, TF, vx_target, vh_target,
        A0_guess_sbagliato, B0_guess_sbagliato,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )

    # Il solver converge (nessuna RuntimeError sollevata sopra) ma trova
    # una soluzione diversa dalla verita' -- il fenomeno delle radici
    # spurie e' riproducibile anche con la fisica corretta.
    assert B_spurio != pytest.approx(B, rel=1e-2)
    # Il segno di B della soluzione spuria e' opposto a quello vero,
    # coerente con la descrizione del fenomeno nella docstring del modulo.
    assert np.sign(B_spurio) != np.sign(B)

    # La soluzione spuria e' comunque una radice valida del residuo (per
    # costruzione fsolve l'ha accettata): lo confermiamo re-integrando.
    riverifica = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A_spurio, B=B_spurio, tf=TF,
        x0=X0_INV, h0=H0_INV, vx0=VX0_INV, vh0=VH0_INV,
    )
    assert riverifica.y[2, -1] == pytest.approx(vx_target, abs=1e-3)
    assert riverifica.y[3, -1] == pytest.approx(vh_target, abs=1e-3)
