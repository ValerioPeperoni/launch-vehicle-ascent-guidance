# Validazione e limiti del progetto

Questo documento raccoglie in un unico posto tutte le approssimazioni
dichiarate in [CLAUDE.md](CLAUDE.md), i loro impatti quantitativi già
misurati durante lo sviluppo, un limite metodologico noto del
risolutore numerico, e un confronto onesto con progetti di riferimento
esterni. Nessun numero qui è una nuova stima: ogni valore è già stato
calcolato e verificato (spesso due volte, in modo indipendente) nel
log dettagliato in [STATUS.md](STATUS.md), a cui si rimanda per il
dettaglio dei calcoli.

## 1. Check di validazione obbligatorio (delta-v)

Il progetto richiede (CLAUDE.md) che il delta-v ideale totale della
traiettoria simulata ricada nel range di plausibilità ~9.1-10.0 km/s
per un'orbita LEO. Con i parametri Falcon 9 usati come caso di
validazione (Step 7) e Isp da vuoto per entrambi gli stadi:

- **Δv ideale totale = 9138.15 m/s** — dentro il range, ma con un
  margine stretto: solo 38 m/s (0.42%) sopra il limite inferiore
  (9100 m/s).
- **Sensibilità con Isp a livello del mare (Isp1 SL) per lo stadio 1
  = 8763.81 m/s** — SOTTO il range. Riportato esplicitamente (non
  omesso) perché la scelta Isp-da-vuoto per lo stadio 1 non è
  arbitraria (lo stadio 1 opera per la maggior parte della salita in
  atmosfera rarefatta, non al livello del mare), ma il check passa con
  un margine che dipende sensibilmente da questa scelta — un lettore
  deve poterlo vedere, non solo il numero che passa.

Questo Δv ideale è un limite di capacità del veicolo (equazione di
Tsiolkovsky, indipendente dalla forma della traiettoria realmente
volata). **Non va confuso** con la domanda, distinta, se lo stato
finale della traiettoria effettivamente integrata corrisponda a
un'orbita stabile: quella domanda è trattata a parte nella Sezione 2
("Niente rotazione terrestre"), dove infatti la traiettoria Falcon 9
del progetto, pur passando questo check di capacità, non raggiunge da
sola un'orbita stabile senza il bonus di rotazione terrestre. I due
fatti non sono in contraddizione: misurano cose diverse (capacità
ideale del veicolo vs. stato cinematico finale raggiunto).

Dettaglio completo del calcolo: STATUS.md, Ciclo 7 (Fase B) e Ciclo 10.

## 2. Approssimazioni dichiarate: impatto misurato

| Approssimazione (da CLAUDE.md) | Impatto misurato | Corretto? |
|---|---|---|
| Atmosfera esponenziale approssimata | Errore ~26% a 11 km, ~2.5% a 25 km rispetto al modello di riferimento (US Standard Atmosphere) | No — dichiarato e accettato fin dall'inizio come semplificazione |
| Cd costante (non varia con Mach) | **Nessun impatto quantificato disponibile.** Il progetto non contiene un confronto Cd-vs-Mach, quindi non si può dire se l'effetto sia piccolo o grande | No — limite non misurato, dichiarato onestamente come tale, non nascosto dietro un "presumibilmente trascurabile" |
| Gravity turn senza sollievo centripeto | 12.2%-17.3% del valore di g alla velocità di fine Stadio 1 (con dati rispettivamente al livello del mare e nel vuoto) | **No** — il termine centripeto è stato aggiunto solo alla fase esoatmosferica (tangente lineare, Ciclo 7 Fase A), non al gravity turn. Asimmetria accettata e dichiarata esplicitamente, non risolta |
| Niente rotazione terrestre | Bonus escluso ≈408.74 m/s (`465.1·cos(28.5°)`, sito equatoriale/subtropicale tipo Cape Canaveral) | Non applicabile (per scope, non un bug) — vedi sotto |
| Resistenza aerodinamica (Cd costante, modello di drag semplificato) | Perdita ideale da drag = 16.7 m/s, contro un range di letteratura generica di ~100-400 m/s | Spiegato, non un errore nascosto — vedi sotto |

### Niente rotazione terrestre: l'esempio più forte del progetto

Questo è il caso in cui "approssimazione dichiarata → impatto misurato
→ conseguenza verificata" si è chiuso nel modo più netto, e vale la
pena raccontarlo per esteso.

Allo Step 9, l'ispezione visiva della traiettoria Falcon 9 validata
allo Step 7 ha mostrato un comportamento sospetto: la quota saliva a
un massimo di ~170 km e poi ridiscendeva a ~148 km a fine bruciamento
(velocità verticale finale negativa). Un calcolo di meccanica orbitale
standard (energia specifica + momento angolare, formule di Curtis) sullo
stato finale di quella traiettoria (h=148.05 km, vx=7755.08 m/s,
vh=-8.24 m/s) ha dato **perigeo ≈ -62.5 km**: sotto la superficie
terrestre, cioè non un'orbita stabile — nonostante la traiettoria
avesse già superato il check di Δv della Sezione 1.

Causa: il risolutore dello Stadio 2 vincolava solo la velocità finale
(vx, vh), lasciando la quota un sottoprodotto non controllato. Corretto
al Ciclo 10 con un nuovo risolutore che vincola esplicitamente quota
finale = target e velocità verticale finale = 0 (volo tangenziale):
il deficit di perigeo si è ridotto a **-22.8 km**. Sono state esplorate
e scartate, con evidenza numerica, altre due leve possibili: un
obiettivo auto-consistente a velocità circolare (migliore risultato
-15.9 km, ancora negativo) e il tempo di bruciamento dello Stadio 2
reso libero come terza incognita (converge sempre al bruciamento
completo, nessun miglioramento — indicando una carenza genuina di
energia/propellente, non un problema di targeting).

Il deficit residuo di ~67 m/s in velocità orizzontale è stato quindi
confrontato con il bonus di rotazione terrestre già escluso per scope
fin dall'inizio del progetto e già quantificato in CLAUDE.md
(408.74 m/s): aggiungendo quel bonus alla velocità orizzontale finale,
il perigeo diventa **esattamente +200.00 km** (il target di quota
stesso, come atteso per costruzione quando la velocità verticale è
già nulla). Il deficit è interamente e quantitativamente coperto dal
bonus, con ampio margine. Verificato due volte in modo indipendente
(orchestratore e critic-ingegnere, con ricalcolo separato).

Non è un difetto di codice o di fisica implementata: è la conseguenza
diretta, ora misurata, di una scelta di scope dichiarata fin
dall'inizio ("niente rotazione terrestre nella prima versione").
Dettaglio completo: STATUS.md, Ciclo 10.

### Drag più basso della letteratura generica: spiegato, non un errore nascosto

La perdita ideale da drag calcolata (16.7 m/s) è molto sotto il range
di letteratura generica citato in CLAUDE.md (100-400 m/s). Non è stato
lasciato come numero isolato sospetto: è stato confrontato con dati
pubblici reali Falcon 9 nella zona di massima pressione dinamica
(Max-Q), che mostrano che il profilo di traiettoria di questo
progetto attraversa quella zona più rapidamente / con un prodotto
Cd·A minore rispetto al veicolo reale completo (fairing, griglie
aerodinamiche, ecc., non modellati). Dettaglio: STATUS.md, Ciclo 9-10.

### Asimmetria non ancora dichiarata altrove: guida a tangente lineare derivata per g costante, integrata con g variabile

Un limite simile nello spirito all'asimmetria del gravity turn (sopra)
esiste anche nella fase esoatmosferica, e va reso esplicito qui. La
legge di guida `tan(θ) = A + B·t` (Step 4) è derivata da Pontryagin
sotto l'ipotesi di gravità linearizzata a valore costante sull'arco di
manovra — è un'approssimazione intrinseca alla tecnica stessa
(Perkins). A partire dal Ciclo 7 (Fase A), però, la dinamica
REALMENTE integrata usa g(h) variabile più il termine di sollievo
centripeto, mentre la forma funzionale della legge di guida (A + B·t)
è rimasta quella derivata per g costante. Dal Ciclo 7 in poi, quindi,
A e B non sono più, in senso stretto, i coefficienti della soluzione
dimostrabilmente ottima di Pontryagin per la fisica effettivamente
simulata: restano una parametrizzazione a due gradi di libertà,
adattata al bersaglio tramite root-finding (Step 6), non una legge
otticamente ottima per l'equazione del moto corretta. Come per il
gravity turn, questa asimmetria guida/dinamica è accettata e non
corretta in questa versione del progetto.

## 3. Limite metodologico: radici multiple nel problema inverso

Il problema inverso della guida a tangente lineare — trovare (A, B)
che portano la traiettoria a un target di velocità terminale dato un
tempo di bruciamento fisso (Step 6) — non è un'approssimazione fisica
ma un limite noto del metodo numerico usato (shooting via
`scipy.optimize.fsolve`/`minimize`).

Il sistema (A, B) → (vx(tf), vh(tf)) ammette almeno due radici
distinte per lo stesso target: una "vera" e una "spuria" con B di
segno opposto (un profilo theta(t) leggermente crescente può produrre
un effetto integrato simile a uno leggermente decrescente su un arco
breve). Esempio numerico concreto (Ciclo 6): con target generato da
A_vero=0.3, B_vero=-0.002, un guess iniziale "ragionevole" ma con
B0=0 converge con residuo a precisione di macchina — cioè sembra una
soluzione corretta — sulla radice **sbagliata** A=0.258, B=+0.002
(segno di B invertito, errore relativo 14% su A, 200% su B). Il
fenomeno è stato confermato sistematico con altre coppie
A_vero/B_vero (0.5/-0.02, 0.4/-0.015), e riscontrato di nuovo, come
previsto, nel contesto più complesso dello staging a metà fase
esoatmosferica (Step 8).

Un guess automatico costruito dalla formula chiusa dello Step 4 (caso
mdot=0) non risolve il problema: converge anch'esso in modo incoerente
sulla radice giusta o sbagliata a seconda del caso, perché il problema
non è "quanto è vicino" il guess in valore assoluto ma **in quale
bacino di attrazione (segno di B)** cade. Mitigazione adottata (non
eliminazione del limite): guess iniziale obbligatorio con
giustificazione fisica esplicita, più un safeguard a doppio guess che
solleva un errore esplicito se due guess ragionevoli ma distinti
convergono a soluzioni diverse, invece di accettare silenziosamente
una delle due.

## 4. Confronto con i progetti di riferimento

Confronto basato su lettura diretta dei README pubblici dei tre
repository citati in CLAUDE.md come riferimento concettuale (nessun
codice o testo copiato, solo parafrasi fattuale). Verificato con due
letture indipendenti (2026-08-18).

| Progetto | Dimensionalità | Stadi | Guida | Validazione con dati reali |
|---|---|---|---|---|
| **Questo progetto** | 2D, piano verticale | Multistadio (con staging anche a metà fase esoatmosferica) | Gravity turn + tangente lineare derivata analiticamente (Perkins/PEG), root-finding per il problema inverso | Sì — Falcon 9, confronto quantitativo con margine esplicito (Sezione 1) |
| [axelstr/gravity_turn_simulation](https://github.com/axelstr/gravity_turn_simulation) | Non dichiarata esplicitamente nel README; verosimilmente piano verticale (nessun riferimento a componenti 3D) | Multistadio (2 stadi, con separazione) | Salita verticale → "programmed turn" ad altitudine data → gravity turn a burn rate costante, con cutoff automatico all'apogeo target + manovre di circolarizzazione | No — nessun confronto numerico con dati reali nel README, solo similarità concettuale dichiarata con Vega/Falcon 9 |
| [bvermeulen/Rocket-and-gravity-turn](https://github.com/bvermeulen/Rocket-and-gravity-turn) | 2D | Monostadio | Controllo ottimo numerico via risolutore CasADi (collocazione diretta), non gravity turn passivo — sofisticazione comparabile alla guida di questo progetto, con metodo diverso | No |
| [b-adkins/pyrocket](https://github.com/b-adkins/pyrocket) | 1D nella versione attuale (0.1.0); 2D è in roadmap (0.2), non implementato | Multistadio (già supportato) | Nessuna guida attiva implementata nella versione attuale (gravity turn è in roadmap) | Non applicabile — scopo dichiaratamente didattico/ricreativo ("intended for aerospace engineering students, amateur rocketeers, and Kerbals") |

Questo confronto **rende più preciso**, senza contraddirlo, il claim
originale di CLAUDE.md secondo cui i progetti individuali reali "usano
tutti punto materiale 2D": resta vero per la dinamica di fondo, ma va
precisato che due dei tre hanno comunque guida attiva sofisticata (uno
multistadio con cutoff automatico e circolarizzazione, l'altro con
controllo ottimo via CasADi) — nessuno dei tre, però, usa una legge di
guida analiticamente derivata e citata in letteratura come questo
progetto (Perkins/PEG), né dichiara un confronto quantitativo con un
benchmark reale con margine esplicito come la Sezione 1 di questo
documento. Gli approcci degli altri due progetti (cutoff euristico su
apogeo, controllo ottimo numerico) sono metodi validi, semplicemente
diversi da quello scelto qui — non è un giudizio di superiorità.

## 5. Stato delle estensioni future dichiarate

Da CLAUDE.md, sezione "Estensioni future dichiarate (fuori scope per
ora)":

- **Rotazione terrestre / 3D completo**: resta fuori scope, come
  dichiarato. Il suo impatto è però ora quantificato e discusso a
  fondo nella Sezione 2 sopra.
- **6-DOF (rotazione del veicolo, momenti aerodinamici)**: resta fuori
  scope, nessuna modifica.
- **Variazione di Cd con il numero di Mach**: resta fuori scope. Come
  notato nella Sezione 2, il suo impatto reale non è quantificato nel
  progetto (nessun dato di confronto disponibile).
- **Ottimizzazione della traiettoria (minimo propellente) oltre alla
  guida a tangente lineare**: resta fuori scope. Nota: lo Step 8
  (staging durante la fase di guida a tangente lineare, con
  ricalcolo di A/B ad ogni segmento) copre un caso — lo staging a
  metà della fase esoatmosferica — che la formulazione originale del
  vincolo "Multistadio" in CLAUDE.md non specificava esplicitamente
  per questa fase; non è però un'ottimizzazione di traiettoria nel
  senso di questa voce (resta guida a target fissato, non minimo
  propellente).
