"""Verifica del modello atmosferico esponenziale dello Step 1.

I valori ISA di riferimento usati come termine di paragone (test 4 e 5)
sono presi da: *US Standard Atmosphere, 1976* (NOAA/NASA/USAF), valori
tabulati anche in Anderson, J.D., *Introduction to Flight*. Sono valori
di un'atmosfera reale/tabulata, usati qui solo per quantificare
l'errore del modello a scala unica del progetto (Curtis Table 8.4,
banda 0-25 km) — non sono la fonte del modello stesso.
"""

import math

import numpy as np
import pytest

from lanciatore.atmosfera import densita
from lanciatore.costanti import RHO0, H_SCALA


def test_densita_a_quota_zero():
    # 1. Valore a quota zero: deve coincidere con RHO0.
    assert densita(0) == pytest.approx(RHO0)


def test_decadimento_monotono():
    # 2. Decadimento monotono con la quota.
    quote = [0, 1000, 5000, 20000, 50000]
    valori = [densita(h) for h in quote]
    for v1, v2 in zip(valori, valori[1:]):
        assert v1 > v2


def test_punto_riferimento_analitico():
    # 3. Punto di riferimento analitico esatto: densita(H_SCALA) == RHO0 / e.
    # Proprieta' matematica dell'esponenziale, indipendente da fonti esterne.
    assert densita(H_SCALA) == pytest.approx(RHO0 / math.e)


def test_confronto_isa_11km():
    # 4. Confronto con ISA reale alla tropopausa (11 km).
    # Fonte del valore di riferimento: US Standard Atmosphere 1976
    # (NOAA/NASA/USAF), rho(11 km) ~= 0.3639 kg/m^3.
    # Il modello a scala singola di questo progetto da' ~= 0.269 kg/m^3,
    # uno scarto di circa il 26%. Questo scarto e' una limitazione nota
    # e accettata del modello a singola scala di altezza (vedi CLAUDE.md
    # e docstring di lanciatore/atmosfera.py): la tolleranza qui e'
    # volutamente ampia (fattore 1.5) solo per verificare che il valore
    # resti nell'ordine di grandezza corretto, NON per nascondere lo
    # scarto ne' per farlo sembrare piu' preciso di quanto sia.
    rho_isa_11km = 0.3639
    rho_modello = densita(11000)
    assert rho_modello == pytest.approx(rho_isa_11km, rel=0.5)
    assert rho_modello / rho_isa_11km < 1.5


def test_confronto_isa_25km():
    # 5. Confronto con ISA reale al bordo della banda del fit (25 km).
    # Fonte del valore di riferimento: US Standard Atmosphere 1976
    # (NOAA/NASA/USAF), rho(25 km) ~= 0.03996 kg/m^3.
    # Qui il modello e' molto piu' preciso (~2.5% di scarto) perche' 25 km
    # e' vicino al bordo della banda su cui Curtis ha fittato H: tolleranza
    # stretta (10%) giustificata per dimostrare la bonta' del fit nella sua
    # banda dichiarata.
    rho_isa_25km = 0.03996
    rho_modello = densita(25000)
    assert rho_modello == pytest.approx(rho_isa_25km, rel=0.10)


def test_asintoto_alta_quota():
    # 6. Ad alta quota la densita' deve restare positiva, piccola, finita.
    rho = densita(100_000)
    assert rho > 0
    assert rho < 1e-4
    assert not math.isnan(rho)
    assert not math.isinf(rho)


def test_quota_negativa_solleva_errore():
    # 7. Quota negativa non fisicamente valida: ValueError.
    with pytest.raises(ValueError):
        densita(-100)


def test_vettorizzazione():
    # 8. Input numpy.ndarray: output ndarray della stessa shape, valori
    # coerenti elemento per elemento rispetto alla chiamata scalare.
    quote = np.array([0, 1000, 5000, 20000, 50000], dtype=float)
    risultati = densita(quote)
    assert isinstance(risultati, np.ndarray)
    assert risultati.shape == quote.shape
    for h, atteso in zip(quote, risultati):
        assert densita(float(h)) == pytest.approx(atteso)
