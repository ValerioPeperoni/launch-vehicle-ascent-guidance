"""Verifica della guida a tangente lineare, fase esoatmosferica (Step 4).

Parametri di test marcati "# PROVVISORIO" (vedi STATUS.md, Ciclo 4,
punto 4): non sono dati di un lanciatore reale (riservato allo Step 6),
ne' una soluzione del problema al contorno A/B-per-target (fuori scope
di questo step, vedi docstring di lanciatore/guida_esoatmosferica.py).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo
modulo (riservato allo Step 6).

IMPORTANTE (criterio 1): i valori attesi in forma chiusa sono
ricalcolati QUI, in codice, con numpy.arcsinh/numpy.sqrt a piena
precisione -- MAI hardcodati come cifre decimali scritte a mano. Il
piano (STATUS.md, Ciclo 4, punto 4 + addendum reviewer punto 8) mostra
che un calcolo a mano di questa stessa formula ha prodotto un errore
reale per cancellazione catastrofica (due asinh vicini moltiplicati per
il fattore grande a0/B=-20000): affidarsi a cifre scritte a mano come
target di pytest.approx sarebbe esattamente l'errore da evitare.
"""

import math

import numpy as np
import pytest

from lanciatore.costanti import G0
from lanciatore.guida_esoatmosferica import (
    derivate_stato_tangente_lineare,
    integra_tangente_lineare,
    risolvi_coefficienti_tangente_lineare,
)

# PROVVISORIO -- Step 4, parametri di test per il confronto in forma
# chiusa (vedi STATUS.md, Ciclo 4, punto 4). mdot=0 e' un caso limite
# scelto apposta: rende a(t) = spinta/m0 ESATTAMENTE costante, quindi
# l'integrale ha forma chiusa esatta (non un'approssimazione).
A = 0.3
B = -0.002  # 1/s
TF = 20.0  # s
M0 = 1000.0  # kg
SPINTA = 40_000.0  # N -> a0 = spinta/m0 = 40 m/s^2
G_COSTANTE = G0  # riuso della costante di progetto (Ciclo 4, punto 4)
X0 = H0 = VX0 = VH0 = 0.0


def test_1_confronto_forma_chiusa():
    # Criterio 1 (il check centrale, vedi STATUS.md Ciclo 4 punto 4):
    # con mdot=0, vx(tf) e vh(tf) calcolati da integra_tangente_lineare
    # devono coincidere con le formule chiuse
    #   vx(tf) - vx0 = (a0/B) * [asinh(A+B*tf) - asinh(A)]
    #   vh(tf) - vh0 = (a0/B) * [sqrt(1+(A+B*tf)^2) - sqrt(1+A^2)] - g*tf
    # entro la tolleranza dell'integratore (rtol=1e-8 usato internamente
    # da integra_tangente_lineare; confronto con rel=1e-6 per margine di
    # accumulo su molti passi, stesso principio gia' giudicato adeguato
    # allo Step 3).
    a0 = SPINTA / M0
    u_f = A + B * TF

    vx_atteso = VX0 + (a0 / B) * (np.arcsinh(u_f) - np.arcsinh(A))
    vh_atteso = (
        VH0
        + (a0 / B) * (np.sqrt(1 + u_f**2) - np.sqrt(1 + A**2))
        - G_COSTANTE * TF
    )

    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    vx_calcolato = risultato.y[2, -1]
    vh_calcolato = risultato.y[3, -1]

    assert vx_calcolato == pytest.approx(vx_atteso, rel=1e-6)
    assert vh_calcolato == pytest.approx(vh_atteso, rel=1e-6)


def test_2_identita_algebrica_diretta():
    # Criterio 2: valutando derivate_stato_tangente_lineare a un
    # istante/stato arbitrario, dvx/dt e dvh/dt restituiti devono
    # coincidere ESATTAMENTE con (spinta/m)*cos(arctan(A+Bt)) e
    # (spinta/m)*sin(arctan(A+Bt))-g (stesso stile del criterio 7 dello
    # Step 3: verifica diretta della formula, indipendente
    # dall'integrazione).
    t_test = 7.3
    y = np.array([123.0, 456.0, 78.0, 9.0, 850.0])
    spinta_test = 55_000.0
    m_test = y[4]

    theta_atteso = np.arctan(A + B * t_test)
    dvx_atteso = (spinta_test / m_test) * np.cos(theta_atteso)
    dvh_atteso = (spinta_test / m_test) * np.sin(theta_atteso) - G_COSTANTE

    derivate = derivate_stato_tangente_lineare(
        t_test, y, A=A, B=B, spinta=spinta_test, mdot=12.0, g_costante=G_COSTANTE
    )

    assert derivate[2] == dvx_atteso
    assert derivate[3] == dvh_atteso


def test_3_monotonia_quota_caso_nominale():
    # Criterio 3: h(t) deve essere strettamente crescente nel caso di
    # test nominale. Plausibilita' fisica gia' verificata nel piano
    # (STATUS.md Ciclo 4, punto 4): dvh/dt parte da
    # 40*sin(16.7 gradi)-9.80665 ~= 1.68 > 0 e resta positivo fino a
    # t=tf, quindi vh(t) > 0 per l'intera durata e h(t) e' monotona
    # crescente. Se questo test fallisse sarebbe un errore nei
    # parametri di test o nell'implementazione, da investigare (non da
    # nascondere), coerente con la regola generale di CLAUDE.md.
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    h = risultato.y[1]
    assert np.all(np.diff(h) > 0)


def test_4_nessun_nan_inf():
    # Criterio 4: nessun NaN/inf su tutta la traiettoria (caso
    # nominale, mdot=0).
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    for componente in risultato.y:
        assert np.all(np.isfinite(componente))


def test_5_massa_costante_quando_mdot_zero():
    # Criterio 5: verifica banale ma diretta che con mdot=0 la massa
    # resti esattamente m0 per tutta la traiettoria (nessun bug che la
    # fa decrescere comunque).
    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=0.0, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    m = risultato.y[4]
    assert m == pytest.approx(M0, rel=0.0, abs=1e-6)


def test_6_caso_con_deplezione_massa_mdot_diverso_zero():
    # Criterio 6 (sanity, non hand-check): con parametri di deplezione
    # realistici (isp=300s coerente con G0, come negli step
    # precedenti), nessun NaN/inf, massa decrescente correttamente. Il
    # confronto analitico esatto resta il caso mdot=0 del criterio 1
    # (la forma chiusa non vale quando la massa depleta davvero, vedi
    # docstring di modulo).
    isp = 300.0  # s
    mdot = SPINTA / (isp * G0)

    # Protezione obbligatoria (addendum reviewer, punto 8 / STATUS.md
    # Ciclo 4 punto 6): verificare ESPLICITAMENTE a monte che mdot*tf
    # resti ben al di sotto di m0, cosi' che spinta/m non rischi mai
    # una singolarita' durante l'integrazione a tf fisso (tf e' fisso e
    # breve in questo scenario di sanity check, nessun evento
    # terminale su massa minima e' presente in questo modulo per
    # costruzione, vedi confine di scope in guida_esoatmosferica.py).
    assert mdot * TF < 0.5 * M0, (
        "mdot*tf non e' sufficientemente piccolo rispetto a m0: rischio di "
        "massa residua troppo bassa entro tf con questi parametri di test."
    )

    risultato = integra_tangente_lineare(
        M0, SPINTA, mdot=mdot, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert risultato.success

    for componente in risultato.y:
        assert np.all(np.isfinite(componente))

    m = risultato.y[4]
    assert np.all(np.diff(m) < 0)  # massa strettamente decrescente
    assert m[-1] == pytest.approx(M0 - mdot * TF, rel=1e-6)


# ---------------------------------------------------------------------------
# Step 6 -- problema inverso: risolvi_coefficienti_tangente_lineare.
#
# Parametri PROVVISORI condivisi da questi test (vedi STATUS.md, Ciclo 6 +
# addendum): mdot REALE (isp=300s, non piu' il caso limite mdot=0 usato per
# la forma chiusa dello Step 4), A_vero/B_vero riusano deliberatamente le
# stesse costanti A, B definite sopra per continuita' col Ciclo 4 (STATUS.md,
# addendum, punto "Valori numerici verificati... per continuita'").
# ---------------------------------------------------------------------------

ISP_STEP6 = 300.0  # s
MDOT_STEP6 = SPINTA / (ISP_STEP6 * G0)  # ~= 13.609 kg/s (vedi addendum)

# Guess raccomandato dall'addendum: A0 = A_vero*0.7, B0 = B_vero*0.7 (stesso
# segno del vero, entro lo stesso bacino di attrazione).
A0_GUESS_CORRETTO = A * 0.7
B0_GUESS_CORRETTO = B * 0.7


def test_7_recupero_non_tautologico_mdot_reale():
    # Criterio non tautologico (STATUS.md Ciclo 6, punto 2.4): il target
    # e' generato integrando AVANTI con A_vero, B_vero e mdot REALE (via
    # integra_tangente_lineare, gia' verificata allo Step 4) -- risultato
    # NUMERICO, non forma chiusa (non esiste con mdot!=0). Il solver deve
    # RECUPERARE A_vero, B_vero partendo da un guess diverso ma nello
    # stesso bacino (A0, B0 scalati 0.7x), non solo azzerare il residuo sul
    # target.
    assert MDOT_STEP6 * TF < 0.5 * M0  # margine di massa ampio (Step 4 criterio 6)

    target = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert target.success
    vx_target = target.y[2, -1]
    vh_target = target.y[3, -1]

    A_risolto, B_risolto = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, MDOT_STEP6, G_COSTANTE, TF, vx_target, vh_target,
        A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )

    # Recupero dei parametri noti (non solo del target), tolleranza stretta.
    assert A_risolto == pytest.approx(A, rel=1e-4)
    assert B_risolto == pytest.approx(B, rel=1e-4)

    # Verifica aggiuntiva: re-integrando con A, B risolti si riottiene il
    # target entro la tolleranza di residuo del solver.
    riverifica = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A_risolto, B=B_risolto,
        g_costante=G_COSTANTE, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert riverifica.y[2, -1] == pytest.approx(vx_target, abs=1e-4)
    assert riverifica.y[3, -1] == pytest.approx(vh_target, abs=1e-4)


def test_8_companion_non_circolare_mdot_zero_formula_chiusa():
    # Test companion (STATUS.md Ciclo 6, addendum, "Test companion non
    # circolare"): il target NON e' generato con integra_tangente_lineare
    # (come nel test 7), ma con la FORMULA CHIUSA ESATTA dello Step 4
    # (asinh/sqrt, stesso stile di test_1_confronto_forma_chiusa sopra),
    # usando mdot=0. Questo rompe la circolarita' generazione-target /
    # funzione-di-ricerca: qui il target non e' nemmeno passato per
    # integra_tangente_lineare prima di essere dato al solver.
    a0 = SPINTA / M0
    u_f = A + B * TF

    vx_target = VX0 + (a0 / B) * (np.arcsinh(u_f) - np.arcsinh(A))
    vh_target = (
        VH0
        + (a0 / B) * (np.sqrt(1 + u_f**2) - np.sqrt(1 + A**2))
        - G_COSTANTE * TF
    )

    A_risolto, B_risolto = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, 0.0, G_COSTANTE, TF, vx_target, vh_target,
        A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )

    assert A_risolto == pytest.approx(A, rel=1e-4)
    assert B_risolto == pytest.approx(B, rel=1e-4)


def test_9_fallimento_esplicito_target_irraggiungibile():
    # Criterio di fallimento (STATUS.md Ciclo 6, addendum): vh_target ben
    # oltre il limite di Tsiolkovsky per questi parametri deve sollevare
    # RuntimeError, non un risultato silenzioso e sbagliato. Il limite si
    # ricalcola qui in codice con math.log (non hardcodato), coerente con
    # lo stile del modulo (nessuna cifra scritta a mano usata come soglia
    # senza essere derivata dal codice).
    m_finale = M0 - MDOT_STEP6 * TF
    dv_ideale_tsiolkovsky = ISP_STEP6 * G0 * math.log(M0 / m_finale)

    vh_target_irraggiungibile = 5000.0  # m/s
    assert vh_target_irraggiungibile > dv_ideale_tsiolkovsky  # verifica esplicita del limite

    vx_target_qualsiasi = 899.267274  # non rilevante: vh da solo gia' impossibile

    with pytest.raises(RuntimeError):
        risolvi_coefficienti_tangente_lineare(
            M0, SPINTA, MDOT_STEP6, G_COSTANTE, TF,
            vx_target_qualsiasi, vh_target_irraggiungibile,
            A0_GUESS_CORRETTO, B0_GUESS_CORRETTO,
            x0=X0, h0=H0, vx0=VX0, vh0=VH0,
        )


def test_10_documentazione_fenomeno_radici_spurie():
    # Documenta il fenomeno scoperto dall'orchestratore (STATUS.md Ciclo 6,
    # addendum): un guess con B0 di segno OPPOSTO al vero converge (residuo
    # piccolo, nessuna RuntimeError) ma su A, B DIVERSI da A_vero, B_vero.
    # Non e' un test che deve "fallire" nel senso pytest -- verifica che il
    # fenomeno e' riproducibile, per proteggere contro una futura modifica
    # che lo nasconda senza che nessuno se ne accorga.
    target = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A, B=B, g_costante=G_COSTANTE, tf=TF,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert target.success
    vx_target = target.y[2, -1]
    vh_target = target.y[3, -1]

    # Guess con B0 di segno OPPOSTO al vero (B_vero = -0.002).
    A0_guess_sbagliato = A * 0.7
    B0_guess_sbagliato = abs(B) * 0.7  # segno invertito rispetto a B_vero

    A_spurio, B_spurio = risolvi_coefficienti_tangente_lineare(
        M0, SPINTA, MDOT_STEP6, G_COSTANTE, TF, vx_target, vh_target,
        A0_guess_sbagliato, B0_guess_sbagliato,
        x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )

    # Il solver converge (nessuna RuntimeError sollevata sopra) ma trova
    # una soluzione diversa dalla verita' -- il fenomeno delle radici
    # spurie e' riproducibile con questi parametri.
    assert B_spurio != pytest.approx(B, rel=1e-2)
    # Il segno di B della soluzione spuria e' opposto a quello vero,
    # coerente con la descrizione del fenomeno nella docstring del modulo.
    assert np.sign(B_spurio) != np.sign(B)

    # La soluzione spuria e' comunque una radice valida del residuo (per
    # costruzione fsolve l'ha accettata): lo confermiamo re-integrando.
    riverifica = integra_tangente_lineare(
        M0, SPINTA, mdot=MDOT_STEP6, A=A_spurio, B=B_spurio,
        g_costante=G_COSTANTE, tf=TF, x0=X0, h0=H0, vx0=VX0, vh0=VH0,
    )
    assert riverifica.y[2, -1] == pytest.approx(vx_target, abs=1e-3)
    assert riverifica.y[3, -1] == pytest.approx(vh_target, abs=1e-3)
