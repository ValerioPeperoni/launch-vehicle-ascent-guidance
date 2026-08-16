"""Verifica del gravity turn, fase atmosferica (Step 3).

Parametri veicolo riusati dallo Step 2 (vedi tests/test_dinamica.py):
non sono dati di un lanciatore reale, quello resta riservato allo
Step 6. Parametri di guida (v_kick, kick_angle_deg) marcati
"# PROVVISORIO": sono scelte di design della legge di guida, non
costanti fisiche di progetto (vedi STATUS.md, Ciclo 3, punto 1.3).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo
modulo (riservato allo Step 6): qui non e' presente ne' multistadio ne'
guida esoatmosferica ne' un caso di validazione con dati reali.
"""

import math

import numpy as np
import pytest

from lanciatore.costanti import CD
from lanciatore.gravita import accelerazione_gravita
from lanciatore.guida import (
    derivate_stato_gravity_turn,
    evento_traiettoria_invalida,
    integra_gravity_turn,
)

# PROVVISORIO — Step 2, parametri veicolo riusati (vedi
# tests/test_dinamica.py per la giustificazione completa di
# plausibilita').
M0 = 50_000.0  # kg
M_VUOTO = 5_000.0  # kg
SPINTA = 800_000.0  # N
ISP = 300.0  # s
DIAMETRO = 3.0  # m
AREA = math.pi * (DIAMETRO / 2) ** 2  # ~= 7.07 m^2

# PROVVISORIO — Step 3, parametri di design della legge di guida (non
# costanti fisiche): vedi STATUS.md, Ciclo 3, punto 1.3.
V_KICK = 50.0  # m/s
KICK_ANGLE_DEG = 2.0  # gradi -> gamma0 = 88 gradi

# Setup comune per test 2 e test 3 (scenario balistico spinta=0, drag=0):
# evita duplicazione del setup di solve_ivp e della scelta di h0/t_span.
Y0_BALISTICO = np.array([0.0, 5000.0, 200.0, np.radians(45.0), M0])
RTOL_BALISTICO = 1e-8
ATOL_BALISTICO = [1e-3, 1e-3, 1e-6, 1e-9, 1e-3]
TSPAN_BALISTICO = (0.0, 20.0)


@pytest.fixture(scope="module")
def risultato_nominale():
    """Integrazione completa (Fase A + Fase B) con i parametri nominali.

    Eseguita una sola volta e riusata dai test che non hanno bisogno di
    condizioni particolari (es. kick esagerato).
    """
    return integra_gravity_turn(
        M0, M_VUOTO, SPINTA, ISP, CD, AREA, V_KICK, KICK_ANGLE_DEG
    )


def test_1_instabilita_equilibrio_gamma_90(risultato_nominale):
    # Criterio 1: a gamma=90 gradi (radianti), dgamma/dt deve essere
    # ~0 per qualunque v>0, indipendentemente da spinta/drag/quota:
    # verifica diretta e analitica (cos(90 gradi)=0) della motivazione
    # del kick angle. Tolleranza a livello di epsilon di macchina (non
    # e' un test approssimato fisicamente, e' zero esatto in aritmetica
    # reale, con un piccolo residuo di floating point per cos(pi/2)).
    gamma_90 = np.radians(90.0)
    for v_test in (10.0, 200.0, 3000.0):
        y = np.array([0.0, 1000.0, v_test, gamma_90, M0])
        derivate = derivate_stato_gravity_turn(
            0.0, y, spinta=SPINTA, isp=ISP, cd=CD, area=AREA, mdot=100.0
        )
        dgamma_dt = derivate[3]
        assert dgamma_dt == pytest.approx(0.0, abs=1e-9)


def test_2_conservazione_velocita_orizzontale_limite_balistico():
    # Criterio 2 (check primario, punto 1.4): con spinta=0, cd=0 (e
    # mdot=0, cosi' la massa non influisce su nulla in questo scenario
    # degenere), v*cos(gamma) deve restare costante entro rtol/atol
    # dell'integratore -- identita' analitica esatta (Culler & Fried),
    # non un margine fisico allargato.
    from scipy.integrate import solve_ivp

    # h0=5000 m e t_span di 20 s (non 60 s): con v0=200 m/s, gamma0=45
    # gradi (vy0~=141 m/s) la traiettoria balistica scenderebbe sotto
    # quota zero entro ~33 s (impatto suolo teorico) -- un intervallo
    # piu' lungo farebbe uscire gli step interni adattivi di RK45 dal
    # dominio h>=0 di accelerazione_gravita/densita (ValueError), che
    # e' un problema di scelta dello scenario di test, non della
    # fisica delle equazioni verificate. h0/t_span scelti con ampio
    # margine (h(t=20s) ~= 3868 m, sempre ben sopra zero).

    risultato = solve_ivp(
        derivate_stato_gravity_turn,
        t_span=TSPAN_BALISTICO,
        y0=Y0_BALISTICO,
        method="RK45",
        args=(0.0, ISP, 0.0, AREA, 0.0),
        rtol=RTOL_BALISTICO,
        atol=ATOL_BALISTICO,
        dense_output=False,
    )
    assert risultato.success

    v = risultato.y[2]
    gamma = risultato.y[3]
    vx = v * np.cos(gamma)

    # Tolleranza coerente con rtol/atol dell'integratore (non fisica
    # allargata): rel=1e-6 e' piu' largo di rtol=1e-8 per assorbire
    # la compounding su una grandezza derivata lungo molti passi, ma
    # resta dell'ordine dell'errore numerico dell'integratore, non un
    # margine di comodo.
    assert vx == pytest.approx(vx[0], rel=1e-6)


def test_3_conservazione_energia_meccanica_specifica():
    # Criterio 3 (check secondario, punto 1.4): stesso scenario del
    # criterio 2, verifica che E = v^2/2 - MU_TERRA/(R_TERRA+h) resti
    # costante entro la stessa tolleranza dell'integratore. Piu' forte
    # del check primario perche' verifica anche la coerenza con la
    # forma specifica di g(h) gia' implementata in gravita.py (non
    # solo una proprieta' generica valida per qualunque g).
    from scipy.integrate import solve_ivp

    from lanciatore.costanti import MU_TERRA, R_TERRA

    # Stesso scenario e stessa motivazione di h0/t_span del criterio 2
    # (vedi commento li'): evita che gli step interni adattivi di RK45
    # escano dal dominio h>=0.

    risultato = solve_ivp(
        derivate_stato_gravity_turn,
        t_span=TSPAN_BALISTICO,
        y0=Y0_BALISTICO,
        method="RK45",
        args=(0.0, ISP, 0.0, AREA, 0.0),
        rtol=RTOL_BALISTICO,
        atol=ATOL_BALISTICO,
        dense_output=False,
    )
    assert risultato.success

    h = risultato.y[1]
    v = risultato.y[2]
    energia = v**2 / 2 - MU_TERRA / (R_TERRA + h)

    assert energia == pytest.approx(energia[0], rel=1e-6)


def test_4_applicazione_kick(risultato_nominale):
    # Criterio 4: subito dopo il kick, lo stato di partenza di Fase B
    # deve avere gamma0 == radians(90 - kick_angle_deg) esatto, x0==0.0
    # esatto, e h0/m0/v0 continui con i valori EFFETTIVI finali della
    # Fase A (non col parametro nominale v_kick, addendum reviewer
    # punto 4) -- test esatto, non approssimato: sono assegnazioni
    # dirette, non il risultato di un'integrazione.
    fase_a = risultato_nominale["fase_verticale"]
    fase_b = risultato_nominale["gravity_turn"]

    h_kick_effettivo, v_kick_effettivo, m_kick_effettivo = fase_a.y[:, -1]
    x0, h0, v0, gamma0, m0 = fase_b.y[:, 0]

    assert x0 == 0.0
    assert gamma0 == np.radians(90.0 - KICK_ANGLE_DEG)
    assert h0 == h_kick_effettivo
    assert v0 == v_kick_effettivo
    assert m0 == m_kick_effettivo


def test_5_evento_traiettoria_invalida_con_kick_esagerato():
    # Criterio 5: con kick_angle_deg=60 (molto oltre il range
    # realistico 1-3 gradi), evento_traiettoria_invalida deve attivarsi
    # (indice 2 nella tupla eventi di Fase B: fine_propellente_2d,
    # impatto_suolo_2d, traiettoria_invalida, velocita_minima).
    # Addendum reviewer punto 2: controllo np.isfinite esteso su TUTTA
    # la traiettoria di entrambe le fasi, non solo nel caso nominale,
    # cosi' un collasso numerico prima dell'attivazione dell'evento non
    # passa inosservato.
    risultato = integra_gravity_turn(
        M0, M_VUOTO, SPINTA, ISP, CD, AREA, V_KICK, kick_angle_deg=60.0
    )
    fase_a = risultato["fase_verticale"]
    fase_b = risultato["gravity_turn"]

    for componente in fase_a.y:
        assert np.all(np.isfinite(componente))
    for componente in fase_b.y:
        assert np.all(np.isfinite(componente))

    t_ev_fine_prop, t_ev_impatto, t_ev_invalida, t_ev_vmin = fase_b.t_events
    assert len(t_ev_invalida) == 1
    assert len(t_ev_fine_prop) == 0


def test_6_caso_nominale_nessun_nan_evento_fine_propellente(risultato_nominale):
    # Criterio 6: con v_kick=50, kick_angle=2 gradi e i parametri
    # veicolo dello Step 2, nessun NaN/inf su nessuna componente di
    # stato di entrambe le fasi; m(t) continua tra le due fasi
    # (verificato gia' nel criterio 4); l'evento che ferma la Fase B
    # deve essere "fine propellente" (indice 0), non uno dei tre eventi
    # difensivi -- se uno di questi scattasse con parametri
    # "ragionevoli" sarebbe un risultato da investigare (bug o
    # combinazione di parametri non fisica), non da nascondere.
    fase_a = risultato_nominale["fase_verticale"]
    fase_b = risultato_nominale["gravity_turn"]

    for componente in fase_a.y:
        assert np.all(np.isfinite(componente))
    for componente in fase_b.y:
        assert np.all(np.isfinite(componente))

    t_ev_fine_prop, t_ev_impatto, t_ev_invalida, t_ev_vmin = fase_b.t_events
    assert fase_b.status == 1
    assert len(t_ev_fine_prop) == 1
    assert len(t_ev_impatto) == 0
    assert len(t_ev_invalida) == 0
    assert len(t_ev_vmin) == 0


def test_7_dgamma_dt_indipendente_da_spinta_e_drag():
    # Criterio 7 (addendum reviewer, punto 3): valutando
    # derivate_stato_gravity_turn con spinta!=0 e drag!=0 (non il caso
    # degenere T=0/D=0) a uno stato arbitrario, la componente
    # dgamma/dt restituita deve coincidere ESATTAMENTE con
    # -(g(h)/v)*cos(gamma), indipendentemente dai valori di
    # spinta/drag passati -- intercetta un eventuale bug che
    # introducesse impropriamente un termine di spinta/drag in
    # dgamma/dt (violerebbe "angolo di attacco nullo" senza che i
    # check di conservazione nel caso degenere se ne accorgano).
    h_test = 5000.0
    v_test = 300.0
    gamma_test = np.radians(37.0)
    y = np.array([1234.0, h_test, v_test, gamma_test, 40_000.0])

    g = accelerazione_gravita(h_test)
    atteso = -(g / v_test) * np.cos(gamma_test)

    for spinta_test, cd_test in ((SPINTA, CD), (SPINTA * 3, CD * 5)):
        derivate = derivate_stato_gravity_turn(
            0.0,
            y,
            spinta=spinta_test,
            isp=ISP,
            cd=cd_test,
            area=AREA,
            mdot=200.0,
        )
        assert derivate[3] == atteso
