"""Modello atmosferico esponenziale approssimato (Step 1).

Formula:

    rho(h) = RHO0 * exp(-h / H_SCALA)      per h >= 0

dove ``h`` e' la quota sul livello del mare in metri.

Fonte del modello: Curtis, H.D., *Orbital Mechanics for Engineering
Students*, Table 8.4 "Exponential Atmospheric Model" — tabella standard
di uso comune in ingegneria aerospaziale che fitta l'atmosfera reale
(US Standard Atmosphere 1976) con un modello esponenziale a tratti
(base altitude, densita' nominale, scala di altezza per banda di
quota). Questo progetto usa la SOLA riga relativa alla banda 0-25 km
(h0 = 0 km, rho0 = 1.225 kg/m^3, H = 7.249 km) come approssimazione
globale a scala unica per tutta l'ascesa, invece del modello a tratti
completo: questa e' la semplificazione dichiarata in CLAUDE.md
("Modello esponenziale approssimato" / "Semplificazione dichiarata
rispetto a US Standard Atmosphere completa").

Fonte alternativa equivalente (nota di confronto, non usata come
valore primario): Anderson, J.D., *Introduction to Flight* —
approssimazione isoterma H = R*T0/g0 con R = 287 J/(kg K),
T0 = 288.15 K, g0 = 9.80665 m/s^2, che da' H ~= 8434 m. Il valore
scelto (7249 m, Curtis) e' preferito perche' e' un fit esplicito
contro dati reali su una banda di quota rilevante per il drag di un
lanciatore, non una singola approssimazione isoterma teorica.

Range di validita' dichiarato: il fit e' accurato principalmente
nella banda 0-25 km. Sopra i ~30 km l'errore relativo cresce, ma la
densita' atmosferica a quelle quote e' comunque cosi' bassa che il
contributo al delta-v da drag e' atteso trascurabile — da verificare
esplicitamente nello Step 6 (validazione), non da assumere qui.
"""

import numpy as np

from lanciatore.costanti import H_SCALA, RHO0


def densita(h):
    """Densita' atmosferica alla quota ``h`` (modello esponenziale a scala unica).

    Parametri
    ---------
    h : float, int o numpy.ndarray
        Quota sul livello del mare, in metri. Deve essere >= 0.
        Input supportati: scalare Python (float/int) o numpy.ndarray.
        Le liste Python native NON sono supportate direttamente: la
        funzione usa operazioni aritmetiche e numpy.exp che si aspettano
        un tipo numerico o un array numpy, non una lista; passare una
        lista solleva TypeError. Convertire esplicitamente con
        numpy.asarray(...) prima della chiamata, se necessario.

    Ritorna
    -------
    float o numpy.ndarray
        Densita' atmosferica in kg/m^3, stessa shape di ``h`` se
        ``h`` e' un array.

    Solleva
    -------
    ValueError
        Se h < 0 (quota non fisicamente valida per un lanciatore
        rispetto al livello del mare; la formula non viene
        estrapolata sotto zero).

    Note
    ----
    Per h molto grande (es. > 100 km) la formula resta valida
    analiticamente (decadimento esponenziale verso 0): non c'e' un
    limite superiore hard-coded ne' un clamp.
    """
    if np.any(np.asarray(h) < 0):
        raise ValueError(
            "Quota h non valida: deve essere >= 0 (livello del mare o superiore), "
            f"ricevuto h={h!r}"
        )
    return RHO0 * np.exp(-h / H_SCALA)
