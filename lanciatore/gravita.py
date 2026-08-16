"""Accelerazione di gravita' variabile con la quota (Step 2).

Formula:

    g(h) = MU_TERRA / (R_TERRA + h)^2

Legge di gravitazione universale per un punto materiale in campo
centrale (Terra sferica non rotante, coerente con l'ipotesi dichiarata
in CLAUDE.md).

Motivazione della scelta g(h) variabile invece di g0 costante (vedi
STATUS.md, Ciclo 2, punto 1.1, per il dettaglio completo):
- MU_TERRA e R_TERRA sono gia' disponibili in costanti.py dallo Step 1
  proprio in vista di questo uso.
- E' l'ipotesi fisicamente piu' corretta e non ha costo implementativo
  aggiuntivo significativo (una divisione in piu' per passo di
  integrazione) rispetto a g0 costante.
- E' coerente con l'uso futuro (Step 4, guida a tangente lineare verso
  un'orbita target: i calcoli di velocita' orbitale ed energia dipendono
  da g(h) preciso). Partire da g0 costante e dover introdurre g(h) in un
  secondo momento creerebbe una discontinuita' concettuale a meta'
  progetto.
- Il check di validazione obbligatorio del delta-v (CLAUDE.md, Step 6)
  richiede la miglior approssimazione disponibile per le perdite
  gravitazionali: usare g0 costante introdurrebbe un errore sistematico
  (a quota orbitale LEO tipica ~200 km, g(h)/g0 ~= 0.94, cioe' ~6% di
  errore) che andrebbe poi isolato e spiegato in fase di validazione
  invece di essere evitato fin da ora.

Nota: g(0) = MU_TERRA / R_TERRA^2 ~= 9.8202 m/s^2, leggermente diverso
da G0 = 9.80665 m/s^2 (scarto relativo ~0.14%). G0 e' un valore
standard convenzionale (include piccoli effetti di
non-sfericita'/rotazione terrestre mediati), mentre g(0) qui e' la pura
gravita' newtoniana di punto materiale a raggio medio: i due valori
sono concettualmente distinti e non ci si aspetta che coincidano
esattamente (vedi tests/test_gravita.py).
"""

import numpy as np

from lanciatore.costanti import MU_TERRA, R_TERRA


def accelerazione_gravita(h):
    """Accelerazione di gravita' alla quota ``h`` (campo centrale, punto materiale).

    Parametri
    ---------
    h : float, int o numpy.ndarray
        Quota sul livello del mare, in metri. Deve essere >= 0.
        Input supportati: scalare Python (float/int) o numpy.ndarray.
        Le liste Python native NON sono supportate direttamente (stesso
        comportamento di lanciatore.atmosfera.densita): convertire
        esplicitamente con numpy.asarray(...) prima della chiamata, se
        necessario.

    Ritorna
    -------
    float o numpy.ndarray
        Accelerazione di gravita' in m/s^2, stessa shape di ``h`` se
        ``h`` e' un array.

    Solleva
    -------
    ValueError
        Se h < 0 (quota non fisicamente valida per un lanciatore
        rispetto al livello del mare), coerentemente con
        lanciatore.atmosfera.densita.
    """
    if np.any(np.asarray(h) < 0):
        raise ValueError(
            "Quota h non valida: deve essere >= 0 (livello del mare o superiore), "
            f"ricevuto h={h!r}"
        )
    return MU_TERRA / (R_TERRA + h) ** 2
