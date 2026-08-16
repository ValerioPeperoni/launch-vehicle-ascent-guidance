"""Verifica dell'accelerazione di gravita' variabile con la quota (Step 2)."""

import math

import numpy as np
import pytest

from lanciatore.costanti import G0, R_TERRA
from lanciatore.gravita import accelerazione_gravita


def test_g0_vicino_a_G0_costante_ma_non_identico():
    # g(0) = MU_TERRA/R_TERRA^2 ~= 9.8202 m/s^2, leggermente diverso da
    # G0 = 9.80665 m/s^2 (scarto relativo atteso ~0.14%). G0 e' un
    # valore standard convenzionale (include piccoli effetti di
    # non-sfericita'/rotazione terrestre mediati), mentre g(0) qui e'
    # la pura gravita' newtoniana di punto materiale a raggio medio:
    # i due valori sono concettualmente distinti, quindi la tolleranza
    # e' stretta ma NON e' un confronto ==.
    g0_calcolato = accelerazione_gravita(0.0)
    assert g0_calcolato == pytest.approx(G0, rel=0.005)
    # Deve comunque essere effettivamente diverso, non un caso in cui
    # per errore accelerazione_gravita ritorna semplicemente G0.
    assert g0_calcolato != G0


def test_decadimento_monotono():
    quote = [0, 1000, 10_000, 100_000, 1_000_000]
    valori = [accelerazione_gravita(h) for h in quote]
    for v1, v2 in zip(valori, valori[1:]):
        assert v1 > v2


def test_mai_negativo_ne_nan():
    quote = np.array([0, 1000, 100_000, 1_000_000, 10_000_000], dtype=float)
    valori = accelerazione_gravita(quote)
    assert np.all(valori > 0)
    assert np.all(np.isfinite(valori))


def test_punto_riferimento_analitico_inverso_quadrato():
    # Proprieta' matematica della legge dell'inverso del quadrato,
    # indipendente da fonti esterne: a h = R_TERRA (distanza dal
    # centro raddoppiata rispetto a h=0), g deve valere circa 1/4 di
    # g(0).
    g_h0 = accelerazione_gravita(0.0)
    g_a_raggio_terra = accelerazione_gravita(R_TERRA)
    assert g_a_raggio_terra == pytest.approx(g_h0 / 4, rel=1e-9)


def test_quota_negativa_solleva_errore():
    with pytest.raises(ValueError):
        accelerazione_gravita(-100)


def test_vettorizzazione():
    quote = np.array([0, 1000, 10_000, 100_000], dtype=float)
    risultati = accelerazione_gravita(quote)
    assert isinstance(risultati, np.ndarray)
    assert risultati.shape == quote.shape
    for h, atteso in zip(quote, risultati):
        assert accelerazione_gravita(float(h)) == pytest.approx(atteso)
