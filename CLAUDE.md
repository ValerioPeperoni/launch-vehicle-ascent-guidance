# Progetto: Simulazione di ascesa e guida di un lanciatore multistadio

## Obiettivo
Simulazione a punto materiale (2D, piano verticale) dell'ascesa propulsa di un
lanciatore multistadio, con guida attiva (non traiettoria scriptata) verso
un'orbita target, con output sia numerico (delta-v, quota/velocita' nel tempo,
budget propellente) sia visivo (animazione della traiettoria).

## Perche' questo scope (non 6-DOF, non 3D con rotazione terrestre)
Verificato con ricerca: i progetti individuali reali con guida attiva verso
un'orbita target (es. axelstr/gravity_turn_simulation, ispirato a Vega/Falcon 9;
bvermeulen/Rocket-and-gravity-turn) usano tutti punto materiale 2D. I progetti
6-DOF completi (RocketPy, CamPyRoS) sono sforzi di team universitari pluriennali
e comunque NON includono guida attiva verso un'orbita — fanno traiettorie
balistiche/non guidate. Lo scope scelto qui e' realistico per un singolo
sviluppatore ed e' anzi piu' specialistico (guida attiva) della maggior parte
dei progetti open source di settore.

## Vincoli tecnici e approssimazioni (non negoziabili senza discussione esplicita)

| Vincolo | Come si implementa | Riferimento |
|---|---|---|
| Dinamica | Punto materiale, 2D (piano verticale), niente rotazione terrestre nella prima versione | Standard per design preliminare di traiettoria |
| Atmosfera | Modello esponenziale approssimato | Semplificazione dichiarata rispetto a US Standard Atmosphere completa |
| Resistenza aerodinamica | Cd costante | Semplificazione dichiarata (in realta' varia con il numero di Mach) |
| Guida, fase atmosferica | Gravity turn (allineamento spinta-velocita' dopo un kick angle iniziale) | Tecnica classica, riferimento Culler et al. 1957 |
| Guida, fase esoatmosferica | Guida a tangente lineare (linear-tangent steering): tan(angolo) = A*t + B | Perkins, "Derivation of Linear-Tangent Steering Laws" — e' l'ottimo matematico esatto per posizione/velocita' dovute alla sola spinta; il PEG (Powered Explicit Guidance, usato dallo Space Shuttle) e' derivato da questo principio |
| Multistadio | Eventi di staging con cambio di massa discontinuo, rilevati come eventi terminali dell'integrazione ODE | Stessa tecnica di gestione eventi gia' usata nel progetto del collasso stellare |
| Caso di validazione | Specifiche pubbliche di un lanciatore reale documentato (es. classe Falcon 9 / Vega) usate come benchmark, non inventate | Stesso principio del catalogo O'Connor & Ott nel progetto precedente |
| **Check di validazione obbligatorio** | Il delta-v totale della traiettoria simulata deve avvicinarsi al range noto ~9.1-10.0 km/s per un'orbita LEO (~7.8 km/s velocita' orbitale + ~1.0-1.5 km/s perdite gravitazionali + ~0.1-0.4 km/s resistenza) | Se lo scarto e' grande, c'e' un errore da isolare, non da nascondere |
| Progetti di riferimento concettuale (NON copiare codice) | axelstr/gravity_turn_simulation, bvermeulen/Rocket-and-gravity-turn, b-adkins/pyrocket | Solo per calibrare scope e approssimazioni ragionevoli |

**Disclaimer esplicito sul check di validazione (aggiunto 2026-08-16,
dopo il primo caso di validazione con dati Falcon 9, Step 7):** il
confronto delta-v vs 9.1-10.0 km/s e' un controllo di **plausibilita' in
un range realistico generico**, NON una riproduzione precisa di uno
specifico lancio reale. Il modello esclude per costruzione il bonus di
velocita' dovuto alla rotazione terrestre (dichiarato fuori scope fin
dall'inizio, vedi riga "Dinamica" sopra) — per un lancio verso est da un
sito equatoriale/subtropicale come Cape Canaveral (28.5°N) questo bonus
vale ≈409 m/s (`465.1 m/s * cos(28.5°)`, velocita' equatoriale di
rotazione proiettata sulla latitudine del sito). Questo valore e' **piu'
grande del margine con cui il check e' stato superato la prima volta**
(38 m/s sopra il limite inferiore del range, Step 7): un confronto
diretto e quantitativo col delta-v di missione di un lancio reale
specifico da un sito equatoriale/subtropicale non sarebbe valido senza
tenerne conto esplicitamente. Il check resta comunque valido come
verifica di ordine di grandezza (il modello non deve produrre un
delta-v totale palesemente sbagliato, es. 5 km/s o 15 km/s), ma non va
presentato ne' interpretato come una riproduzione accurata delle
prestazioni di un lanciatore reale specifico.

## Estensioni future dichiarate (fuori scope per ora, non implementare senza richiesta esplicita)
- Rotazione terrestre / 3D completo
- 6-DOF (rotazione del veicolo, momenti aerodinamici)
- Variazione di Cd con il numero di Mach
- Ottimizzazione della traiettoria (minimo propellente) oltre alla guida a tangente lineare

## Sub-agenti del progetto (in .claude/agents/)
- planner: scompone il progetto in step, aggiorna STATUS.md
- reviewer: secondo controllo critico su piani e report
- coder: implementa il codice dello step corrente
- critic-ingegnere: verifica il codice contro la tabella dei vincoli sopra, incluso il check del delta-v
- optimizer: pulizia codice, stile, ridondanze (modello economico)
- reporter: scrive il resoconto giornaliero sintetico (modello economico)

## Regole di lavoro
- Prima di implementare un nuovo step, usa la plan mode per proporre l'approccio.
- Filtra sempre l'output numerico verboso (log di simulazione) prima che arrivi
  nel contesto principale.
- Ogni ciclo giornaliero: leggi STATUS.md, esegui il prossimo step, aggiorna
  STATUS.md, genera il report con il subagent reporter, fermati e attendi
  conferma dell'utente.
- Se un test di convergenza/verifica numerica da' un risultato inatteso, NON
  aggiustare tolleranze per farlo passare senza prima capire la causa fisica —
  stesso principio gia' seguito nel progetto del collasso stellare (vedi lo
  Step 4 di quel progetto come esempio di comportamento corretto).
- Nessun push pubblico su GitHub senza conferma esplicita dell'utente.
