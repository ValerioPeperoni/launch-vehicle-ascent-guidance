"""Verifica delle costanti fisiche dello Step 1.

Valori da letteratura, non calcolati: i test verificano tipo numerico
ed ordine di grandezza, non un valore esatto arbitrario (tranne dove
il valore stesso e' una definizione, es. G0).
"""

from lanciatore.costanti import G0, R_TERRA, MU_TERRA, RHO0, H_SCALA


def test_g0():
    assert isinstance(G0, float)
    assert 9.7 < G0 < 9.9


def test_r_terra():
    assert isinstance(R_TERRA, float)
    # Raggio medio terrestre atteso attorno a 6371 km.
    assert 6_300_000 < R_TERRA < 6_400_000


def test_mu_terra():
    assert isinstance(MU_TERRA, float)
    assert 3.9e14 < MU_TERRA < 4.0e14


def test_rho0():
    assert isinstance(RHO0, float)
    # Densita' aria a livello del mare, ISA: attorno a 1.225 kg/m^3.
    assert 1.1 < RHO0 < 1.3


def test_h_scala():
    assert isinstance(H_SCALA, float)
    # Scala di altezza Curtis Table 8.4 (banda 0-25 km): 7.249 km.
    assert 7000 < H_SCALA < 7500
