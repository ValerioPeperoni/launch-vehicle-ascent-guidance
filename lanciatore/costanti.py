"""Costanti fisiche del progetto (Step 1).

Ogni costante e' documentata con valore, unita' SI e fonte. Le costanti
qui presenti sono quelle necessarie a partire dallo Step 1; altre
costanti (Cd, parametri di stadio) verranno aggiunte negli step
successivi quando servono, per evitare dead code.
"""

# Gravita' standard terrestre.
# Valore: 9.80665 m/s^2. Definizione internazionale (3a CGPM, 1901).
G0 = 9.80665

# Raggio medio terrestre.
# Valore: 6 371 000 m. Fonte: IUGG. Coerente con l'assunzione di Terra
# sferica non rotante dichiarata in CLAUDE.md.
R_TERRA = 6_371_000.0

# Parametro gravitazionale standard terrestre (GM).
# Valore: 3.986004418e14 m^3/s^2. Valore WGS84, riportato anche in
# Curtis, "Orbital Mechanics for Engineering Students" e in Vallado,
# "Fundamentals of Astrodynamics and Applications".
MU_TERRA = 3.986004418e14

# Densita' atmosferica a livello del mare (ISA standard).
# Valore: 1.225 kg/m^3. Coincide con l'entry a quota 0 km di Curtis
# Table 8.4 ("Exponential Atmospheric Model").
RHO0 = 1.225

# Scala di altezza del modello atmosferico esponenziale a scala unica.
# Valore: 7249 m (7.249 km). Fonte: Curtis, "Orbital Mechanics for
# Engineering Students", Table 8.4, riga relativa alla banda di quota
# base h0 = 0 km. Vedi lanciatore/atmosfera.py per il dettaglio del
# modello e per una fonte alternativa (isoterma, Anderson) citata come
# nota di confronto.
H_SCALA = 7249.0

# Coefficiente di resistenza aerodinamica (Cd), costante per l'intera
# traiettoria (semplificazione dichiarata in CLAUDE.md: "Resistenza
# aerodinamica: Cd costante" rispetto alla reale dipendenza dal numero
# di Mach, esplicitamente fuori scope).
# Valore: 0.3. Rappresentativo per un lanciatore slanciato (fusoliera
# cilindrica con ogiva), a meta' del range 0.2-0.5 indicato in
# CLAUDE.md. Coerente con valori tipici per corpi affusolati in regime
# transonico/supersonico riportati in Anderson, "Introduction to
# Flight", e in Sutton, "Rocket Propulsion Elements".
CD = 0.3
