# Stato del progetto — aggiornato ad ogni ciclo

## Step
- [x] Step 1: Setup ambiente, modello atmosferico esponenziale, costanti fisiche
- [x] Step 2: Dinamica a punto materiale (spinta, gravita', drag) — singolo stadio, ascesa verticale pura, senza guida (caso di test piu' semplice)
- [x] Step 3: Gravity turn (fase atmosferica)
- [x] Step 4: Guida a tangente lineare (fase esoatmosferica) — verifica contro derivazione analitica nota
- [x] Step 5: Multistadio, fase atmosferica (eventi di staging, cambio massa discontinuo — scope deciso 2026-08-16: solo verticale/gravity turn, non tangente lineare, vedi Step 7 sotto)
- [x] Step 6: Guida a tangente lineare — problema inverso (root-finding
  A/B per un target di velocita' terminale, deplezione di massa reale)
  — **inserito 2026-08-16**, scorporato da quello che era lo Step 6
  originale (vedi nota sotto), perche' la validazione con dati reali non
  puo' produrre una traiettoria completa senza prima saper risolvere A/B
  per un target. Root-finding via shooting (scipy.optimize), nessuna
  forma chiusa disponibile con mdot reale (a differenza del caso
  mdot=0 gia' verificato allo Step 4).
- [x] Step 7: Caso di validazione con dati reali di un lanciatore pubblico + confronto delta-v vs benchmark ~9.1-10.0 km/s
  - **Nota dell'utente (2026-08-16, ORA RISOLTA — vedi nuovo Step 6
    sopra):** lo Step 4 ha implementato la tangente lineare con A/B DATI
    (integra, non risolve). Questo step (ex-Step 6) richiede il problema
    INVERSO, ora scorporato nel nuovo Step 6 dedicato — qui si riusa il
    risolutore del nuovo Step 6, non lo si reimplementa.
- [x] Step 8: Estensione — staging durante la fase di guida a tangente
  lineare (aggiunto 2026-08-16, NON opzionale per il progetto finale,
  solo rimandato nell'ordine — vedi decisione utente nel Ciclo 5).
  Root-finding di A/B per segmento (si appoggia al problema inverso, ora
  Step 6) + verifica di continuita' di stato tra uno stadio e il
  successivo quando lo staging avviene a meta' della fase esoatmosferica
  (i coefficienti A/B del segmento in corso non sono piu' validi dopo un
  cambio discontinuo di massa/spinta, vanno ricalcolati per il nuovo
  segmento).
  - **Nota dell'utente (2026-08-16, da NON riscoprire da capo):** lo
    Step 6 ha trovato che il problema inverso A/B ha radici multiple/
    spurie per un singolo segmento (vedi Ciclo 6, addendum). Con più
    segmenti concatenati con continuita' tra loro (staging a metà della
    fase esoatmosferica), è molto probabile che lo stesso fenomeno si
    ripresenti, probabilmente amplificato (più incognite, più bacini di
    attrazione possibili). Il safeguard costruito allo Step 6 (due guess
    vicini che devono concordare, altrimenti RuntimeError esplicito) è
    probabilmente riusabile qui, non da reinventare — verificarlo
    esplicitamente prima di progettare un meccanismo nuovo da zero.
- [x] Step 9: Visualizzazione (traiettoria numerica + animazione)
  - **Scoperta 2026-08-17, RISOLTA al Ciclo 10 (2026-08-18):**
    ispezionando visivamente i grafici prodotti da questo step, emerso
    che la quota nel Segmento 2 dello Step 7 sale a un massimo di
    ~170km poi RIDISCENDE a ~148km a fine bruciamento (vh diventa
    negativo). Approfondito dal critic-ingegnere e verificato
    indipendentemente dall'orchestratore (leggi orbitali standard,
    nessuna nuova fisica): lo stato finale della traiettoria Step 7
    (h=148.05km, vx=7755.08, vh=-8.24 m/s) NON corrisponde a un'orbita
    LEO stabile — perigeo calcolato ≈-62.5 km (sotto la superficie
    terrestre), apogeo ≈148.3 km. **Risolto al Ciclo 10:** risolutore
    corretto per vincolare quota+verticalità (non più solo velocità
    orizzontale), scarto ridotto a -22.8km; esplorate e scartate con
    evidenza numerica altre due leve (obiettivo a velocità circolare
    auto-consistente: -15.9km, ancora negativo; tempo di bruciamento
    libero come terza incognita: converge sempre al bruciamento
    completo, nessun miglioramento); **chiuso definitivamente
    verificando che il deficit residuo (~67 m/s) è interamente spiegato
    dal bonus di rotazione terrestre (409 m/s) già escluso per scope fin
    dall'inizio — con quel bonus il perigeo diventa +200km esatto.**
    Verificato due volte in modo indipendente (orchestratore +
    critic-ingegnere). Dettaglio completo: STATUS.md Ciclo 10,
    disclaimer aggiornato in CLAUDE.md (2026-08-18).
- [x] Step 10: Validazione e documentazione dei limiti (VALIDATION.md, confronto concettuale con i progetti di riferimento)
- [ ] Step 11: Pulizia, documentazione, README, preparazione per GitHub

## Log cicli
(ogni ciclo aggiunge una riga qui: data, step completato, note)

---

### 2026-08-15 — Ciclo 1 (planner): piano dettagliato Step 1

**Step:** 1 — Setup ambiente, modello atmosferico esponenziale, costanti fisiche
**Stato:** pianificazione completata. Contesto gia' fatto (non da ripianificare): venv creato,
requirements.txt (numpy/scipy/matplotlib/pytest), repo git locale, .gitignore, agenti in .claude/agents/.
**Prossimo:** passare al sub-agente coder per l'implementazione secondo questo piano.

#### 1. Struttura dei file da creare

```
lanciatore/
├── __init__.py         # pacchetto vuoto (eventualmente __version__ = "0.1.0")
├── costanti.py          # tutte le costanti fisiche del progetto, con unita' SI e fonte in docstring/commento per ciascuna
└── atmosfera.py          # funzione densita(h) col modello esponenziale approssimato

tests/
├── test_costanti.py      # verifica valori/ordini di grandezza delle costanti
└── test_atmosfera.py     # verifica del modello atmosferico (criteri al punto 4)

conftest.py                # (root del progetto, vuoto) — garantisce che pytest trovi il pacchetto
                            # lanciatore/ in layout flat senza bisogno di installarlo in modo editable
```

Nota per il coder: NON creare `tests/__init__.py`, pytest con rootdir al livello del progetto
individua i test senza bisogno di renderla una sub-package; se in fase di esecuzione pytest non trova
`lanciatore`, aggiungere `conftest.py` vuoto in root (forza pytest ad aggiungere la root a sys.path)
prima di ricorrere a soluzioni piu' invasive (pyproject.toml con package config, pip install -e).

#### 2. Costanti fisiche (`lanciatore/costanti.py`)

Costanti necessarie ORA per lo Step 1:

| Nome variabile | Simbolo | Valore | Unita' | Fonte |
|---|---|---|---|---|
| `G0` | g0 | 9.80665 | m/s^2 | Gravita' standard, definizione internazionale (3a CGPM, 1901) |
| `R_TERRA` | R_T | 6 371 000 | m | Raggio medio terrestre (IUGG), coerente con l'assunzione di Terra sferica non rotante dichiarata in CLAUDE.md |
| `MU_TERRA` | mu | 3.986004418e14 | m^3/s^2 | Parametro gravitazionale standard terrestre GM, valore WGS84 (stesso valore riportato in Curtis, *Orbital Mechanics for Engineering Students*, e in Vallado, *Fundamentals of Astrodynamics*) |
| `RHO0` | rho0 | 1.225 | kg/m^3 | Densita' atmosferica a livello del mare, ISA standard — coincide con l'entry a quota 0 km di Curtis Table 8.4 (vedi punto 3) |
| `H_SCALA` | H | 7249 | m (7.249 km) | Scala di altezza — Curtis Table 8.4, riga a quota base h0 = 0 km (vedi punto 3 per dettagli) |

Costanti PREVISTE per step successivi (da NON implementare ora, solo per consapevolezza — se il coder
le anticipa per errore, il critic-ingegnere deve segnalarlo come scope creep):
- `CD` — coefficiente di resistenza aerodinamica costante (Step 2), valore tipico atteso 0.2–0.5 per un lanciatore slanciato, da fissare col caso di validazione dello Step 6.
- Costanti di ogni stadio (massa a vuoto, massa propellente, spinta, Isp) — Step 5/6, dipendono dal
  lanciatore di riferimento scelto per la validazione.

Non includere `R` (costante specifica dei gas per l'aria, 287 J/(kg·K)) e `T0` (288.15 K, temperatura
ISA al livello del mare) come costanti attive nel modulo: servono solo come riferimento per *derivare*
H (vedi punto 3), non sono usate direttamente nella dinamica. Il coder puo' citarle nella docstring di
`H_SCALA` come nota di derivazione alternativa, ma non aggiungerle come variabili module-level inutilizzate
(evitare dead code che l'optimizer dovrebbe poi rimuovere).

#### 3. Modello atmosferico esponenziale (`lanciatore/atmosfera.py`)

**Formula:**

```
rho(h) = RHO0 * exp(-h / H_SCALA)      per h >= 0
```

dove `h` e' la quota sul livello del mare in metri.

**Fonte esplicita (da citare nella docstring del modulo/funzione):**
Curtis, H.D., *Orbital Mechanics for Engineering Students*, Table 8.4 "Exponential Atmospheric Model" —
tabella standard di uso comune in ingegneria aerospaziale che fitta l'atmosfera reale (US Standard
Atmosphere 1976) con un modello esponenziale a tratti (base altitude, densita' nominale, scala di
altezza per banda di quota). Il progetto usa la SOLA riga relativa alla banda 0–25 km
(h0 = 0 km, rho0 = 1.225 kg/m^3, H = 7.249 km) come approssimazione globale a scala unica per tutta
l'ascesa, invece del modello a tratti completo — questa e' esattamente la semplificazione dichiarata
in CLAUDE.md ("Modello esponenziale approssimato" / "Semplificazione dichiarata rispetto a US Standard
Atmosphere completa").

Fonte alternativa equivalente (da menzionare come nota, non da usare come valore primario, per
trasparenza sulla scelta): Anderson, J.D., *Introduction to Flight* — approssimazione isoterma
H = R*T0/g0 con R = 287 J/(kg K), T0 = 288.15 K, g0 = 9.80665 m/s^2, che da' H ≈ 8434 m. Il valore
scelto (7249 m, Curtis) e' preferito perche' e' un fit esplicito contro dati reali su una banda di
quota rilevante per il drag di un lanciatore, non una singola approssimazione isoterma teorica.

**Range di validita' dichiarato:** il fit e' accurato principalmente nella banda 0–25 km. Sopra i
~30 km l'errore relativo cresce, ma la densita' atmosferica a quelle quote e' comunque cosi' bassa che
il contributo al delta-v da drag e' trascurabile — da verificare pero' esplicitamente nello Step 6
(non assumerlo, controllarlo nel check di validazione delta-v).

**Comportamento ai bordi (da specificare nella funzione, decisione di design per il coder):**
- `h < 0`: solleva `ValueError` (quota non fisicamente valida per un lanciatore rispetto al livello
  del mare; non estrapolare la formula sotto zero).
- `h` molto grande (es. >100 km): la formula resta valida analiticamente (decadimento esponenziale
  verso 0), nessun clamp necessario — non c'e' un limite superiore hard-coded.
- La funzione deve accettare sia uno scalare Python/float sia un array numpy per `h` (verra' usata
  vettorizzata dentro l'integratore ODE nello Step 2+): usare `numpy.exp`, non `math.exp`.

#### 4. Criteri di verifica (`tests/test_atmosfera.py`)

1. **Valore a quota zero:** `densita(0) == RHO0` (1.225 kg/m^3), confronto esatto o con tolleranza
   numerica minima (es. `pytest.approx`).
2. **Decadimento monotono:** per una sequenza di quote crescenti (es. 0, 1000, 5000, 20000, 50000 m),
   la densita' deve essere strettamente decrescente.
3. **Punto di riferimento analitico esatto:** `densita(H_SCALA) == RHO0 / e` (proprieta' matematica
   dell'esponenziale, indipendente da qualunque fonte esterna — buon test di correttezza
   dell'implementazione in se').
4. **Confronto con ISA reale a 11 km (tropopausa):** valore ISA reale ≈ 0.3639 kg/m^3; il modello
   a scala singola da' ≈ 0.269 kg/m^3 (~26% di scarto). Il test deve verificare che il valore sia
   nell'ordine di grandezza corretto con tolleranza ampia (es. entro un fattore 1.5), E il test deve
   contenere un commento esplicito che documenta questo scarto come limitazione nota del modello a
   scala singola — NON restringere la tolleranza per farla sembrare piu' precisa di quanto sia, e NON
   allargarla oltre il necessario per nascondere il problema. Questo e' il caso concreto menzionato in
   CLAUDE.md: uno scarto atteso e capito va documentato, non nascosto.
5. **Confronto con ISA reale a 25 km (bordo banda del fit):** valore ISA reale ≈ 0.03996 kg/m^3; il
   modello da' ≈ 0.0389 kg/m^3 (~2.5% di scarto, molto piu' preciso qui perche' e' vicino al bordo
   della banda su cui Curtis ha fittato H). Tolleranza stretta qui (es. entro 10%) e' giustificata e
   serve a dimostrare che il fit e' buono nella sua banda dichiarata.
6. **Asintoto ad alta quota:** a h = 100 000 m la densita' deve essere un numero positivo molto
   piccolo (es. < 1e-4 kg/m^3), mai negativa, mai NaN/inf.
7. **Quota negativa:** `densita(-100)` deve sollevare `ValueError` (`pytest.raises`).
8. **Vettorizzazione:** passare un `numpy.ndarray` di quote e verificare che l'output sia un ndarray
   della stessa shape con valori coerenti elemento per elemento rispetto alla chiamata scalare.

**`tests/test_costanti.py`:** verifica che ciascuna delle 5 costanti dello Step 1 esista, abbia il tipo
numerico atteso, e sia nell'ordine di grandezza corretto (es. `G0` tra 9.7 e 9.9, `MU_TERRA` tra 3.9e14
e 4.0e14, ecc.) — non serve altro per questo step, sono valori da letteratura, non calcolati.

#### 5. Task granulari per il coder (ordine consigliato)

1. Creare `lanciatore/__init__.py` (vuoto).
2. Creare `lanciatore/costanti.py` con le 5 costanti della tabella al punto 2, ciascuna con commento
   inline che riporta valore/unita'/fonte.
3. Creare `lanciatore/atmosfera.py` con la funzione `densita(h)` secondo la formula e le regole ai
   bordi del punto 3, docstring con la fonte (Curtis Table 8.4) citata per esteso.
4. Creare `conftest.py` vuoto in root, se necessario per la scoperta del pacchetto da parte di pytest.
5. Creare `tests/test_costanti.py` e `tests/test_atmosfera.py` con i criteri del punto 4.
6. Eseguire `pytest` e riportare solo il riepilogo (pass/fail count), non il log verboso — filtrare
   prima che l'output arrivi al contesto principale, come da regola di lavoro in CLAUDE.md.

#### 6. Note per il critic-ingegnere

- Verificare che la fonte di H (Curtis Table 8.4) sia effettivamente citata per iscritto nella
  docstring di `atmosfera.py`, non solo qui in STATUS.md.
- Verificare che il test a 11 km NON sia stato reso artificialmente stretto o largo per farlo passare
  senza commento esplicativo — deve essere leggibile perche' lo scarto e' accettato.
- Segnalare esplicitamente se ritiene che lo scarto del 26% a 11 km meriti di essere rivalutato in uno
  step successivo (es. modello a due scale di altezza) — per ora e' una scelta consapevole, ma va
  ridiscussa se lo Step 6 (validazione delta-v) mostra un contributo di drag anomalo.
- Questo step non tocca ancora Cd, massa, spinta: nessun check di delta-v e' pertinente qui.

#### 7. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **approvabile con modifiche minori**. Modifiche recepite:

1. **Fonte dei valori ISA di riferimento nei test 4-5 (mancava):** i valori
   rho(11km) ~= 0.3639 kg/m^3 e rho(25km) ~= 0.03996 kg/m^3 usati come termine
   di paragone vanno citati in `test_atmosfera.py` con fonte esplicita:
   *US Standard Atmosphere, 1976* (NOAA/NASA/USAF), valori tabulati anche in
   Anderson, *Introduction to Flight* — non solo la fonte del modello va
   citata, anche quella dei benchmark di confronto.
2. **Criterio "niente dead code" chiarito:** G0, R_TERRA, MU_TERRA restano in
   `costanti.py` anche se non consumate da `atmosfera.py` in questo step,
   perche' sono costanti di progetto (servono agli step 2/5/6, e sono
   verificate ora da `test_costanti.py`). R e T0 restano invece SOLO in
   docstring perche' sono costanti di derivazione intermedia, non costanti di
   progetto: non verranno mai usate direttamente altrove nel codice.
3. **Input non numpy/scalare:** `densita(h)` NON supporta liste Python native
   (solleverebbe TypeError con `numpy.exp` su una lista) — da scrivere
   esplicitamente nella docstring: input supportati float, int, numpy.ndarray.
4. **`conftest.py`:** crearlo SEMPRE (non condizionale) come deliverable
   definito, vuoto, in root — costo zero, evita un giro di debug sulla
   scoperta del pacchetto da parte di pytest.
5. **`__init__.py`:** vuoto, senza `__version__` per ora (verra' introdotto
   allo Step 9 in fase di preparazione per GitHub/packaging).

#### 8. Esito ciclo (coder + critic-ingegnere + optimizer)

- **coder**: implementati `lanciatore/__init__.py`, `lanciatore/costanti.py`,
  `lanciatore/atmosfera.py`, `conftest.py`, `tests/test_costanti.py`,
  `tests/test_atmosfera.py`, secondo il piano + addendum. Nessuna ambiguita'
  residua trovata.
- **critic-ingegnere**: pytest rieseguito in modo indipendente (13 passed).
  Vincolo "Atmosfera: modello esponenziale approssimato" **verificato**,
  fonte citata nel codice (non solo in STATUS.md). Tutti gli altri vincoli
  della tabella CLAUDE.md sono "non pertinenti a questo step" (nessuna
  dinamica/guida/staging/delta-v implementata ancora) — nessuna violazione.
  Nota aperta per i cicli futuri: il check delta-v vs 9.1-10.0 km/s va
  eseguito esplicitamente quando disponibile (atteso Step 6), non omesso.
- **optimizer**: un solo intervento cosmetico (range(len()) -> zip() in
  test_decadimento_monotono). 13 test confermati verdi dopo la modifica.

**Step 1: COMPLETATO.** Prossimo step proposto: Step 2 (dinamica a punto
materiale, singolo stadio, ascesa verticale pura senza guida).

---

### 2026-08-15 — Ciclo 2 (planner): piano dettagliato Step 2

**Step:** 2 — Dinamica a punto materiale (spinta, gravita', drag) — singolo
stadio, ascesa verticale pura, senza guida (caso di test piu' semplice).
**Stato:** pianificazione completata, in attesa di passaggio al coder.

Questo step introduce la dinamica vera e propria ma nella forma piu' semplice
possibile: moto puramente verticale (1 grado di liberta' spaziale: quota h),
un solo stadio, nessuna guida (l'angolo spinta-velocita' e' fisso a 90°
rispetto all'orizzonte per definizione, non c'e' ancora ne' gravity turn ne'
tangente lineare — quelli sono gli Step 3/4). Serve come base numerica
verificabile analiticamente (Tsiolkovsky) prima di aggiungere la complessita'
della guida.

#### 1. Punti aperti risolti

**1.1 Gravita': costante o variabile con la quota?**

**Decisione: variabile con la quota**, `g(h) = MU_TERRA / (R_TERRA + h)^2`
(legge di gravitazione universale per un punto materiale, campo centrale).

Motivazione:
- `MU_TERRA` e `R_TERRA` sono gia' in `costanti.py` dallo Step 1 proprio in
  vista di questo uso (vedi nota nel log del Ciclo 1: "servono agli step
  2/5/6").
- E' l'ipotesi fisicamente piu' corretta e non ha alcun costo implementativo
  aggiuntivo significativo (una divisione in piu' per passo di integrazione).
- E' coerente con l'uso futuro (Step 4, guida a tangente lineare verso
  un'orbita: la velocita' orbitale target e i calcoli di energia dipendono da
  g(h) preciso). Iniziare con g0 costante e poi dover passare a g(h) in uno
  step successivo introdurrebbe una discontinuita' concettuale da rigiustificare
  a meta' progetto.
- Il check di validazione obbligatorio (delta-v vs 9.1-10.0 km/s, Step 6)
  richiede che le perdite gravitazionali siano calcolate con la miglior
  approssimazione disponibile con lo stesso costo implementativo di quella
  costante: usare g0 costante introdurrebbe un errore sistematico (g0 e'
  valida solo al livello del mare; a quota orbitale LEO tipica ~200 km,
  g(h)/g0 ~= 0.94, quindi ~6% di errore sistematico sulle perdite
  gravitazionali) che poi andrebbe isolato e spiegato in fase di validazione
  invece di essere evitato da subito.

Sanity check numerico incluso nei test (vedi punto 1.4): `g(0)` calcolato con
la formula deve essere vicino a `G0` ma NON identico — `g(0) = MU_TERRA /
R_TERRA^2 ≈ 9.8202 m/s^2` contro `G0 = 9.80665 m/s^2` (scarto relativo
~0.14%). Questo scarto e' atteso e noto: `G0` e' un valore standard
convenzionale (include piccoli effetti di non-sfericita'/rotazione terrestre
mediati), mentre `g(0)` qui e' la pura gravita' newtoniana di punto materiale
a raggio medio. Il test deve verificare che i due valori coincidano entro una
tolleranza stretta ma non nulla (es. 0.5%) — NON deve essere un test
`==` esatto, e il commento nel test deve spiegare la ragione dello scarto (non
lasciarlo come "magic tolerance" non spiegata).

**1.2 Parametri del veicolo per questo caso di test**

**Decisione: (b) numeri di test plausibili e dichiaratamente provvisori**,
non ancora dati di un lanciatore reale. Il caso di validazione con dati
pubblici reali e' esplicitamente riservato allo Step 6 (vedi CLAUDE.md,
tabella vincoli, riga "Caso di validazione"). Usare gia' ora dati reali di
un lanciatore specifico creerebbe due problemi: (a) anticiperebbe lo Step 6
senza che siano ancora implementati staging, guida e drag realistico, quindi
il confronto sarebbe comunque privo di senso a questo stadio; (b) se poi allo
Step 6 servisse aggiustare qualcosa nella struttura dati per adattarla al
lanciatore di validazione scelto, si rischierebbe di dover ritoccare un
riferimento "reale" gia' pubblicato, confondendo cosa e' benchmark e cosa e'
placeholder.

I parametri di test NON vanno in `costanti.py` (che contiene solo costanti
fisiche di progetto, sourced, valide per tutto il progetto) ma restano
locali a `tests/test_dinamica.py`, con commento esplicito
`# PROVVISORIO — Step 2, sostituito da dati reali allo Step 6` sopra la loro
definizione. `lanciatore/dinamica.py` resta agnostico rispetto al veicolo: la
funzione di integrazione riceve massa iniziale, massa a vuoto, spinta, Isp,
Cd e area come parametri, non li hard-codea.

Valori di test proposti e giustificazione di plausibilita' (ordini di
grandezza per un singolo stadio a propellente liquido kerolox, confrontabili
con un primo stadio orbitale reale ma NON copiati da uno specifico):

| Parametro | Valore | Plausibilita' |
|---|---|---|
| Massa iniziale `m0` | 50 000 kg | Ordine di grandezza di un primo stadio orbitale piccolo/medio (per confronto, un vero primo stadio Falcon 9 e' ~500 t: qui si sceglie deliberatamente un valore piu' piccolo e "tondo" per un caso di test semplice e veloce da integrare, non per imitare un lanciatore specifico) |
| Massa a vuoto `m_vuoto` | 5 000 kg | Rapporto di massa m0/m_vuoto = 10 → frazione di propellente 90%, nell'intervallo tipico (85-92%) per uno stadio a propellenti liquidi ad alte prestazioni |
| Massa propellente `m_prop` | 45 000 kg | = m0 - m_vuoto |
| Spinta `T` | 800 000 N | Da T/W0 = T/(m0*G0) = 800000/490332 ≈ 1.63: rapporto spinta/peso al lift-off nel range tipico 1.2-2.0 per un primo stadio |
| Impulso specifico `Isp` | 300 s | Tipico per un motore kerolox a livello del mare (RP-1/LOx, range realistico ~280-310 s) |

Grandezze derivate (da riportare in docstring/commento nel test, per
trasparenza, non da hard-codare come costanti separate):
- Portata massica (costante durante il bruciamento): `mdot = T/(Isp*G0) ≈
  271.9 kg/s`
- Tempo di bruciamento atteso: `t_burn = m_prop/mdot ≈ 165.5 s` (ordine di
  grandezza plausibile per un primo stadio, tipicamente 150-200 s)
- Delta-v ideale di Tsiolkovsky (nessuna gravita', nessun drag): `Isp*G0*ln(m0/m_vuoto)
  ≈ 6774 m/s` — usato come limite superiore nel check del punto 1.5.

**Nota importante per il critic-ingegnere:** questi valori NON vanno
confrontati con il benchmark delta-v 9.1-10.0 km/s del CLAUDE.md — quel
check e' esplicitamente riservato allo Step 6 con dati reali e con lo stadio
completo (incluse guida e eventualmente multistadio). In questo step al
massimo ci si aspetta un delta-v raggiunto nettamente inferiore (ascesa
verticale pura = perdite gravitazionali massime possibili, nessuna orbita
raggiunta, un solo stadio).

**1.3 Cd e area di riferimento**

Nuova costante in `costanti.py`:

| Nome | Valore | Fonte/motivazione |
|---|---|---|
| `CD` | 0.3 | Valore rappresentativo per un lanciatore slanciato (fusoliera cilindrica con ogiva), a meta' del range 0.2-0.5 indicato in CLAUDE.md. Coerente con valori tipici riportati per corpi affusolati in regime transonico/supersonico in Anderson, *Introduction to Flight*, e Sutton, *Rocket Propulsion Elements*. Costante per l'intera traiettoria per costruzione (CLAUDE.md: "Cd costante" e' una semplificazione dichiarata rispetto alla reale dipendenza dal numero di Mach — NON implementare la variazione con Mach, e' esplicitamente fuori scope) |

L'area di riferimento (sezione trasversale del veicolo), a differenza di
`CD`, e' una proprieta' specifica del veicolo, non una costante di
modellazione valida per l'intero progetto: resta quindi tra i parametri di
test provvisori in `tests/test_dinamica.py`, non in `costanti.py`, per lo
stesso principio del punto 1.2. Valore di test proposto: diametro 3.0 m →
`area = pi * (3.0/2)^2 ≈ 7.07 m^2` (ordine di grandezza plausibile per un
primo stadio orbitale reale, es. Falcon 9 ~3.7 m, Vega ~3.0 m — usato solo
come riferimento di scala, non come dato del lanciatore di validazione).

**1.4 Metodo di integrazione e criterio di arresto**

- **Vettore di stato** (caso 1D verticale): `y = [h, v, m]` (quota, velocita'
  verticale, massa corrente). Non serve ancora una componente orizzontale:
  l'ascesa e' puramente verticale in questo step per costruzione.
- **Equazioni del moto:**
  - `dh/dt = v`
  - `dv/dt = (T - D - m*g(h)) / m`, con `D = 0.5 * rho(h) * CD * area * v^2 *
    sign(v)` (drag sempre opposto al moto; per v >= 0 in questo step si
    riduce a `D = 0.5*rho(h)*CD*area*v^2` con verso negativo — la forma con
    `sign(v)` e' scritta comunque in modo generale per riusabilita' futura,
    quando in Step 3+ la velocita' potra' avere componenti non puramente
    verticali)
  - `dm/dt = -mdot`, con `mdot = T/(Isp*G0)` costante durante il bruciamento
    (nota: la portata massica e' definita per convenzione internazionale
    usando `G0` costante anche se la gravita' del moto usa `g(h)` variabile —
    sono due grandezze fisicamente distinte, non vanno confuse)
- **Integratore:** `scipy.integrate.solve_ivp`, metodo esplicito di default
  (RK45) — il sistema non e' rigido (stiff) per questo caso, coerente con
  l'uso di RK45/RK4 nei progetti di riferimento concettuale citati in
  CLAUDE.md. Tolleranze: `rtol=1e-8`, `atol` come vettore componente per
  componente (non uno scalare unico) perche' h (~10^3-10^5 m), v (~10^2-10^3
  m/s) e m (~10^3-10^4 kg) hanno scale molto diverse — es.
  `atol=[1e-3, 1e-6, 1e-3]`.
- **Evento terminale "fine propellente":** funzione evento
  `evento_fine_propellente(t, y) = y[2] - m_vuoto` (massa corrente meno massa
  a vuoto), con `terminal=True`, `direction=-1` (si annulla scendendo da
  positivo a zero). Questo e' esattamente il pattern richiesto da CLAUDE.md
  per gli eventi di staging (Step 5) — qui viene gia' usato in forma
  "degenere" con un solo stadio, cosi' la struttura del codice non cambia
  quando arrivera' lo staging vero, cambia solo il numero di eventi/stadi.
  **Nessuna fase di volo balistico (coasting) dopo il burnout va simulata in
  questo step**: l'integrazione si ferma esattamente all'evento di fine
  propellente, per costruzione (il caso di test e' definito come "ascesa
  propulsa pura").
- **Evento di sicurezza "impatto suolo"** (`h(t) = 0` con `direction=-1`,
  `terminal=True`): da includere gia' ora come evento aggiuntivo difensivo,
  anche se con i parametri scelti (T/W0 ≈ 1.63) non dovrebbe mai attivarsi —
  serve a stabilire fin da ora il pattern "piu' eventi terminali coesistono"
  che sara' necessario allo Step 5 con lo staging multiplo, e protegge da
  integrazioni patologiche se in futuro qualcuno passa parametri non validi.
- **Controllo di validita' iniziale:** la funzione di integrazione deve
  sollevare **`ValueError`** (non `assert`, che e' disattivabile con
  `python -O` e quindi inadatto a validare parametri forniti da
  utente/test) se `T <= m0 * g(0)` prima di integrare — altrimenti il
  "lift-off" non decolla per definizione fisica e il test di monotonia del
  punto 1.5 non avrebbe senso.
- **Caso limite a t=0 per l'evento "impatto suolo":** con `h0=0.0` il
  trigger dell'evento (`h(t)=0`) coincide col valore iniziale. Non deve
  scattare spuriamente: `solve_ivp` rileva un evento sul cambio di segno
  della funzione evento *durante* l'integrazione, non sul valore isolato a
  t=0, e con `v0=0`, `dv/dt(0)>0` la quota si allontana da zero
  immediatamente dopo l'istante iniziale. Da verificare esplicitamente nel
  test (punto 1.5, criterio 5) controllando che `t_events` dell'evento
  "impatto_suolo" sia vuoto.
- **Struttura di file suggerita per la funzione:**
  `integra_ascesa_verticale(m0, m_vuoto, spinta, isp, cd, area, h0=0.0,
  v0=0.0, t_max=...)` in `lanciatore/dinamica.py`, che ritorna l'oggetto
  risultato di `solve_ivp` (o i suoi array `t`, `y` estratti) — non
  nascondere il risultato grezzo, il test deve poter accedere a
  `t_events`/`status` per verificare quale evento ha fermato l'integrazione
  (vedi punto 1.5).

**1.5 Criteri di verifica per questo step**

Nessun confronto con il benchmark delta-v reale (Step 6) e nessuna guida
(Step 3/4) in questo step. Controlli di sanita' fisica da implementare in
`tests/test_dinamica.py`:

1. **Massa:** `m(0) == m0` esatto; `m(t)` strettamente decrescente per tutta
   la durata del bruciamento; valore finale `m(t_fine) == m_vuoto` (entro
   tolleranza numerica dell'integratore, es. `pytest.approx` con tolleranza
   coerente con `atol` usato).
2. **Monotonia quota/velocita':** `h(t)` strettamente crescente per tutta la
   traiettoria (conseguenza diretta di `v > 0` dopo il lift-off). `v(t)`
   strettamente crescente finche' la spinta netta E' POSITIVA, cioe' finche'
   `T > D(t) + m(t) * g(h(t))` — **correzione rispetto alla bozza iniziale**:
   la condizione deve includere il drag `D(t)`, non solo il peso, perche'
   `dv/dt = (T - D - m*g(h))/m` (punto 1.4): usare `T > m*g` da solo
   testerebbe una condizione necessaria ma non sufficiente e potrebbe far
   fallire il test per un picco di drag anche in assenza di bug. Per
   evitare di duplicare in modo disallineato l'implementazione del drag nel
   test, il criterio va derivato direttamente dal segno di `dv/dt` calcolato
   con le stesse equazioni/array che l'integratore produce (non da un
   ricalcolo indipendente di `D(t)` nel test). Con questa correzione, un
   eventuale fallimento del test e' inequivocabilmente un bug (nel drag o
   nei parametri), non un caso limite da valutare caso per caso — quindi
   **non** va comunque "aggiustato" allargando tolleranze, coerente con la
   regola generale di CLAUDE.md su risultati inattesi.
3. **Nessun NaN/inf:** `np.isfinite(...)` su tutti gli array `h`, `v`, `m`
   dell'intera soluzione.
4. **Confronto con Tsiolkovsky ideale (limite superiore, non test stretto):**
   `v_finale < Isp * G0 * ln(m0/m_vuoto)` (≈ 6774 m/s con i parametri di
   test) — deve essere strettamente minore (per via delle perdite
   gravitazionali e di drag), ma il test verifica solo la direzione della
   disuguaglianza, non un valore target preciso. Riportare anche nel test
   (come commento o print in caso di fallimento) il rapporto
   `v_finale/v_ideale`: essendo un'ascesa verticale pura senza gravity turn,
   ci si aspetta un rapporto relativamente basso (perdite gravitazionali
   massime possibili, `cos(0)=1` per l'intera durata) — questo e' un
   comportamento atteso di questo caso di test degenere, non un difetto da
   correggere.
5. **Evento di terminazione corretto:** verificare tramite `result.status` /
   `result.t_events` che l'integrazione si sia fermata per l'evento "fine
   propellente" E, esplicitamente, che `t_events` dell'evento "impatto
   suolo" sia **vuoto** (non dedurlo solo indirettamente dal fatto che
   "fine_propellente" e' quello popolato) — con questi parametri l'evento
   di impatto non deve mai attivarsi, ne' a t=0 (vedi nota caso limite al
   punto 1.4) ne' durante l'ascesa; se si attiva e' un bug o un parametro di
   test errato, da investigare.
6. **Condizioni iniziali:** `h(0) == 0`, `v(0) == 0` esatti.

`tests/test_gravita.py` (nuovo modulo, analogo a `test_atmosfera.py`):
verifica `g(0)` vicino a `G0` entro 0.5% (vedi punto 1.1, con commento che
spiega lo scarto atteso), decadimento monotono di `g(h)` per quote crescenti,
mai negativo/NaN, e un punto di riferimento analitico (es. a `h = R_TERRA`,
`g` deve essere circa un quarto di `g(0)`, proprieta' dell'inverso del
quadrato — buon test di correttezza dell'implementazione indipendente da
fonti esterne).

#### 2. Struttura dei file da creare

```
lanciatore/
├── gravita.py            # nuova: accelerazione_gravita(h) = MU_TERRA/(R_TERRA+h)^2
└── dinamica.py            # nuova: equazioni del moto (ascesa verticale pura,
                            #   niente guida) + integra_ascesa_verticale(...)
                            #   con solve_ivp ed eventi terminali (fine
                            #   propellente + impatto suolo di sicurezza)

tests/
├── test_gravita.py         # verifica g(h), criteri del punto 1.5 (ultimo blocco)
└── test_dinamica.py         # parametri veicolo di test PROVVISORI (punto 1.2/1.3)
                              #   + criteri di verifica del punto 1.5 (1-6)
```

`costanti.py`: aggiunta della sola costante `CD` (punto 1.3). Nessun'altra
modifica ai file dello Step 1.

#### 3. Task granulari per il coder (ordine consigliato)

1. Aggiungere `CD = 0.3` a `lanciatore/costanti.py` con commento
   valore/fonte come da tabella al punto 1.3 (NON aggiungere l'area di
   riferimento qui, resta nel test).
2. Creare `lanciatore/gravita.py` con `accelerazione_gravita(h)` secondo la
   formula del punto 1.1, docstring con motivazione della scelta
   gravita'-variabile-vs-costante (riassunta dal punto 1.1), stesso stile di
   `atmosfera.py` (supporto scalare + numpy.ndarray, `ValueError` per h < 0
   coerentemente con `atmosfera.densita`).
3. Creare `lanciatore/dinamica.py`:
   a. funzione delle derivate di stato (equazioni del punto 1.4);
   b. funzione evento "fine propellente" e funzione evento "impatto suolo";
   c. funzione `integra_ascesa_verticale(...)` che assembla tutto e chiama
      `solve_ivp` con gli eventi, tolleranze e controllo iniziale
      T > m0*g(0) del punto 1.4.
4. Creare `tests/test_gravita.py` con i criteri del punto 1.5 (blocco
   `test_gravita.py`).
5. Creare `tests/test_dinamica.py`: definire i parametri di test provvisori
   del punto 1.2/1.3 con il commento `# PROVVISORIO` richiesto, poi i 6
   criteri di verifica del punto 1.5.
6. Eseguire `pytest`, riportare solo il riepilogo pass/fail (non il log
   numerico verboso), come da regola di lavoro CLAUDE.md.

#### 4. Note per il critic-ingegnere

- Verificare che la scelta g(h) variabile sia effettivamente implementata
  (non g0 costante "per ora") e che la motivazione sia nel codice
  (docstring di `gravita.py`), non solo qui in STATUS.md.
- Verificare che i parametri veicolo in `test_dinamica.py` siano
  esplicitamente marcati come provvisori con un commento visibile nel
  codice, e che NON siano stati inseriti (per errore o comodita') in
  `costanti.py` insieme a `CD`.
- Vincolo CLAUDE.md pertinente a questo step: "Resistenza aerodinamica: Cd
  costante" — verificare che `CD` sia effettivamente usato come costante
  singola (nessuna dipendenza da Mach/velocita' introdotta).
- Vincolo "Multistadio: eventi di staging come eventi terminali
  dell'integrazione ODE" — non ancora pertinente allo staging vero (Step 5),
  ma verificare che il pattern con `solve_ivp` + evento terminale "fine
  propellente" sia gia' impostato in modo generalizzabile (stessa firma di
  funzione evento riusabile), come richiesto esplicitamente in questo piano.
- **Nessun check di delta-v vs benchmark 9.1-10.0 km/s e' pertinente in
  questo step** (riservato allo Step 6) — non segnalarlo come mancante.
- Punto di attenzione esplicito (vedi punto 1.5.2): se il test di monotonia
  di `v(t)` rivela un calo locale per drag/max-Q con questi parametri di
  test, e' un risultato da investigare (causa fisica: parametri poco
  realistici? errore nel segno del drag?), non da nascondere allargando la
  tolleranza — coerente con la regola di lavoro CLAUDE.md su risultati
  inattesi nei test di verifica numerica.

#### 5. Esito ciclo (coder + critic-ingegnere + optimizer)

- **coder**: implementati `lanciatore/gravita.py` (accelerazione_gravita),
  `lanciatore/dinamica.py` (derivate_stato, eventi fine_propellente/
  impatto_suolo, integra_ascesa_verticale), `CD=0.3` in costanti.py,
  `tests/test_gravita.py`, `tests/test_dinamica.py` (parametri veicolo
  marcati `# PROVVISORIO`). 25/25 test verdi (13 Step1 + 12 Step2). Nessuna
  ambiguita' residua.
- **critic-ingegnere**: pytest e metriche ricalcolate in modo indipendente
  (coincidono col report del coder). Vincoli "punto materiale 2D",
  "atmosfera esponenziale", "Cd costante" **verificati**. Pattern eventi
  ODE per staging **verificato come preparazione corretta** (non ancora
  testabile in pieno con un solo stadio). Parametri veicolo confermati
  provvisori e assenti da costanti.py. Delta-v parziale di questo step
  (v_finale=5090 m/s, rapporto v/v_ideale=0.7514) confrontato esplicitamente
  col benchmark 9.1-10.0 km/s: scarto atteso e spiegato (ascesa verticale
  pura senza guida = perdite gravitazionali massime, singolo stadio,
  parametri provvisori) — NON e' il check di validazione vero, che resta
  riservato allo Step 6. Nessuna violazione.
- **optimizer**: nessuna correzione necessaria, codice gia' pulito. 25 test
  confermati verdi.

**Step 2: COMPLETATO.** Prossimo step proposto: Step 3 (gravity turn, fase
atmosferica).

---

### 2026-08-15 — Ciclo 3 (planner): piano dettagliato Step 3

**Step:** 3 — Gravity turn (fase atmosferica).
**Stato:** pianificazione completata, in attesa di passaggio al coder.

Questo step introduce la prima forma di guida attiva: dopo un breve tratto di
ascesa verticale pura (fisicamente identica allo Step 2, riusata cosi' com'e'),
si applica un kick angle iniziale e si passa a un moto 2D nel piano verticale
in cui la spinta e' sempre allineata alla velocita' (angolo di attacco nullo)
e l'assetto evolve secondo l'equazione classica del gravity turn. Restiamo per
costruzione in fase atmosferica per l'intera durata di questo step: nessuna
logica di passaggio alla guida a tangente lineare (Step 4) viene introdotta
qui.

#### 1. Punti aperti risolti

**1.1 Formulazione dello stato 2D — polare/flight-path per la fase di gravity
turn.**

**Decisione:** per la fase di gravity turn (Fase B, vedi punto 1.3) lo stato
e' `y = [x, h, v, gamma, m]` (distanza a terra, quota, modulo velocita',
angolo di rotta rispetto all'orizzonte locale — flight-path angle —, massa).
`gamma` e' portato internamente in **radianti** per tutta l'integrazione
(conversione da gradi solo nei parametri di ingresso leggibili da umani, es.
il kick angle).

Motivazione (richiesta esplicitamente: relazione con solve_ivp e col
riferimento analitico):
- L'equazione di guida del gravity turn (punto 1.2, Culler & Fried 1957) e'
  scritta *direttamente* in termini di `dgamma/dt`. Portare `gamma` come
  variabile di stato esplicita significa che l'integratore calcola
  esattamente la grandezza che va verificata contro la formula di
  riferimento — nessuna ricostruzione indiretta di `gamma` (es. via
  `atan2(vh, vx)`) e nessun problema di ramo/quadrante nell'arcotangente.
  Questo rende il confronto col riferimento analitico (punto 1.4) diretto e
  trasparente per il critic-ingegnere: si puo' ricalcolare
  `-g(h)/v * cos(gamma)` sui vettori di stato prodotti da `solve_ivp` e
  confrontarlo punto per punto con la derivata numerica di `gamma(t)`.
- Una formulazione cartesiana (`x, h, vx, vh`) non evita comunque il calcolo
  di `gamma`: per proiettare la spinta lungo la velocita' (allineamento
  spinta-velocita' richiesto dal gravity turn) serve comunque `gamma =
  atan2(vh, vx)` a ogni passo, quindi la formulazione cartesiana
  aggiungerebbe un calcolo trigonometrico extra senza rimuovere alcuna
  criticita' numerica in questo caso specifico.
- Il potenziale vantaggio della forma cartesiana (evitare singolarita' a
  `v=0` o gestire `gamma` che attraversa lo zero/negativo, es. in fasi
  balistiche dopo l'apogeo) non e' rilevante per lo scope di questo step:
  la Fase B parte gia' con `v = v_kick > 0` (mai da `v=0`, vedi punto 1.3) e
  resta in ascesa propulsa atmosferica, senza attraversare l'apogeo o
  ricadere. Se in step futuri (fuori scope qui) servisse gestire coasting
  balistico con inversione di `gamma`, la forma cartesiana potrebbe essere
  rivalutata in quel momento — non anticiparla ora.
- Coerente con lo stile dei progetti di riferimento concettuale citati in
  CLAUDE.md (es. axelstr/gravity_turn_simulation, che integra in forma
  `v, gamma, h, x`, secondo la tecnica classica citata anche da Culler &
  Fried) — solo per calibrare lo scope, non copiato.

La **Fase A** (ascesa verticale pre-kick, vedi punto 1.3) resta invece nello
stato 1D `[h, v, m]` gia' definito ed *effettivamente riusato* da
`lanciatore/dinamica.py` (Step 2) — nessuna riformulazione necessaria li',
vedi punto 1.6.

**1.2 Equazione di guida del gravity turn.**

**Equazione dichiarata esplicitamente (da citare nella docstring del
codice, non solo qui):**

```
dv/dt     = (T - D) / m - g(h) * sin(gamma)
dgamma/dt = -(g(h) / v) * cos(gamma)
dh/dt     = v * sin(gamma)
dx/dt     = v * cos(gamma)
dm/dt     = -mdot
```

con `T` spinta sempre allineata alla velocita' (angolo di attacco nullo per
ipotesi di gravity turn — nessuna componente di portanza/normale), `D =
0.5 * rho(h) * CD * area * v^2` (drag lungo la velocita', verso opposto,
`v >= 0` per costruzione in Fase B), `g(h)` come gia' implementato in
`lanciatore/gravita.py` (riusato, non ridefinito — vedi punto 1.6), `x`
distanza orizzontale nel piano verticale locale (approssimazione a piano
piatto non rotante, coerente con l'ipotesi 2D gia' dichiarata in CLAUDE.md;
la curvatura terrestre nel legame tra `x` e la geometria del campo di
gravita' resta fuori scope, come gia' per lo Step 2).

**Fonte esplicita (obbligatoria, da citare nel codice):** Culler, G.J. e
Fried, B.D., "General Theory of Gravity Turn Trajectories", *Journal of
Applied Physics* / *ARS Journal*, 1957 — riferimento gia' indicato in
CLAUDE.md per questa tecnica di guida. E' la fonte dell'equazione
`dgamma/dt = -(g/v)*cos(gamma)`, derivata imponendo che l'angolo di attacco
sia nullo (spinta e velocita' sempre allineate) durante tutta la fase
atmosferica propulsa.

**1.3 Kick angle iniziale.**

Il gravity turn puro ha un equilibrio (instabile) esattamente a `gamma =
90°`: a quell'angolo `cos(gamma) = 0`, quindi `dgamma/dt = 0` per costruzione
— senza una perturbazione iniziale il veicolo continuerebbe a salire in
verticale indefinitamente nel modello (fatto verificabile direttamente: vedi
test dedicato al punto 4). Serve quindi un kick angle applicato una tantum,
tecnica standard citata anche da Culler & Fried e usata in tutti i progetti
di riferimento concettuale di CLAUDE.md.

**Decisione — sequenza a due fasi:**
- **Fase A (ascesa verticale pura, pre-kick):** identica per costruzione
  allo Step 2 (stato `[h, v, m]`, `gamma` implicitamente 90° per l'intera
  durata, `dx/dt = 0`). Riusa **direttamente** `derivate_stato` di
  `lanciatore/dinamica.py` (vedi punto 1.6), integrata con `solve_ivp` fino
  a un nuovo evento terminale "raggiunta velocita' di kick"
  (`y[1] - v_kick`, `terminal=True`, `direction=+1`), definito localmente
  nel nuovo modulo (punto 1.6) — NON nel corpo di `dinamica.py`, per non
  toccare codice gia' testato dallo Step 2.
- **Kick (istantaneo, a fine Fase A):** al raggiungimento di `v = v_kick` si
  ridefinisce lo stato per la Fase B come `x0 = 0` (la Fase A e' puramente
  verticale, quindi la distanza a terra e' zero fino a quel momento),
  `h0 = h_kick`, `v0 = v_kick` (continuita' del modulo della velocita'),
  `gamma0 = 90° - kick_angle_deg` (discontinuita' voluta e unica variabile
  che cambia istantaneamente), `m0 = m_kick`. L'idealizzazione "istantanea"
  e' un'approssimazione standard in letteratura (il beccheggio fisico reale
  avviene su una finestra temporale breve ma non nulla, a bassa pressione
  dinamica appena dopo il decollo — qui semplificata a un salto discreto,
  coerente con l'approccio dei progetti di riferimento citati in CLAUDE.md).
- **Valori di test proposti (parametri di guida, non costanti fisiche —
  vedi motivazione strutturale al punto 1.6):** `v_kick = 50 m/s`,
  `kick_angle_deg = 2.0°` (→ `gamma0 = 88°`). Plausibilita' numerica coi
  parametri veicolo provvisori dello Step 2 (T/W0 ≈ 1.63, accelerazione
  netta iniziale ≈ 0.63*g ≈ 6.2 m/s^2): `v_kick = 50 m/s` viene raggiunto in
  circa 8 s, a una quota di circa 200 m — ordine di grandezza plausibile per
  un pitch-over di un lanciatore reale (basso, a pressione dinamica ancora
  contenuta). L'angolo di 2° e' nell'intervallo tipico (1-3°) citato nei
  progetti di riferimento concettuale per questo tipo di simulazione.

**Nota strutturale importante:** `v_kick` e `kick_angle_deg` sono parametri
del **design della legge di guida**, non costanti fisiche di progetto e non
proprieta' del veicolo in senso stretto — restano quindi parametri
espliciti delle funzioni del nuovo modulo (nessun default hard-coded dentro
`lanciatore/guida.py`), con i valori numerici di test marcati
`# PROVVISORIO` in `tests/test_guida.py`, esattamente come i parametri
veicolo dello Step 2 (stesso principio, stessa motivazione: sono scelte di
progettazione/ottimizzazione specifiche del singolo lanciatore, il modulo di
guida deve restare agnostico).

**1.4 Fonte analitica per la verifica (obbligatoria).**

Confermato l'approccio proposto, con un secondo check complementare che lo
rafforza:

**Check primario — conservazione della velocita' orizzontale nel limite
balistico in vuoto (T=0, D=0):** con spinta e drag nulli, dalle due equazioni
del punto 1.2:

```
d(v*cos(gamma))/dt = cos(gamma)*dv/dt - v*sin(gamma)*dgamma/dt
                    = cos(gamma)*(-g*sin(gamma)) - v*sin(gamma)*(-(g/v)*cos(gamma))
                    = -g*sin(gamma)*cos(gamma) + g*sin(gamma)*cos(gamma)
                    = 0
```

quindi `v*cos(gamma)` (= `vx`) e' **esattamente conservato**, analiticamente,
per qualunque legge `g(h)` (costante o variabile con la quota: `g` si
elide algebricamente, la proprieta' non dipende dalla forma di `g(h)`).
Derivazione ottenuta direttamente dalle due equazioni gia' citate al punto
1.2 (Culler & Fried) — nessuna fonte esterna aggiuntiva necessaria, e'
matematica pura applicata a un'equazione gia' sourced.

**Check secondario (aggiunto in questo piano, cross-check indipendente sulla
stessa Fase B, stesso limite T=0/D=0):** conservazione dell'energia
meccanica specifica `E = v^2/2 - MU_TERRA/(R_TERRA + h)`, coerente con il
modello di gravita' effettivamente implementato in `gravita.py`
(`g(h) = MU_TERRA/(R_TERRA+h)^2 = d/dh[MU_TERRA/(R_TERRA+h)]`, quindi
`g(h)` e' esattamente il gradiente del potenziale usato in `E`). Fonte:
Curtis, *Orbital Mechanics for Engineering Students*, definizione
dell'energia meccanica specifica per un campo centrale (gia' testo di
riferimento del progetto, citato per `MU_TERRA` fin dallo Step 1). Questo
secondo check e' piu' forte del primo perche' verifica anche la coerenza
tra l'equazione di gravity turn e la *forma specifica* di `g(h)` gia'
implementata e testata (non solo una proprieta' generica valida per
qualunque `g`), quindi intercetta anche eventuali errori di segno/scala
nell'uso di `accelerazione_gravita(h)` dentro le nuove equazioni.

Entrambi i check vanno eseguiti impostando `spinta=0` e `cd=0` (o `area=0`)
nella funzione delle derivate della Fase B, integrando su un arco di alcuni
secondi/minuti con condizioni iniziali arbitrarie ma fisicamente sensate
(es. `v=200 m/s`, `gamma=45°`, `h=1000 m`), e verificando che
`v*cos(gamma)` ed `E` restino costanti entro la tolleranza numerica
dell'integratore (non entro una tolleranza fisica ampia: sono identita'
esatte, quindi la tolleranza accettabile e' quella di `rtol`/`atol` scelti
per `solve_ivp`, non un margine arbitrario).

**1.5 Confine di questo step.**

Confermato: lo Step 3 integra la Fase B (gravity turn) fino all'evento
"fine propellente" dello stadio corrente (stesso pattern dello Step 2, m_vuoto
riusato), oppure fino a un evento difensivo (impatto suolo o traiettoria
invalida, vedi punto 1.6) — **non** viene introdotta alcuna condizione di
uscita basata sulla quota (es. "quota di Karman ~100 km" o simili) come
proxy per il passaggio alla guida a tangente lineare: quella decisione
(quando e come cambiare legge di guida) appartiene esplicitamente allo
Step 4 e non va anticipata qui, nemmeno in forma di placeholder o TODO nel
codice.

**1.6 Riuso vs duplicazione con dinamica.py.**

**Decisione: nuovo modulo `lanciatore/guida.py`**, che:
- **riusa senza modifiche** `lanciatore.atmosfera.densita` e
  `lanciatore.gravita.accelerazione_gravita` (import diretto, nessuna
  reimplementazione della fisica di drag/gravita');
- **riusa senza modifiche** `lanciatore.dinamica.derivate_stato` per la
  Fase A (import diretto della funzione, passata cosi' com'e' a `solve_ivp`
  con un evento terminale nuovo definito localmente in `guida.py` — vedi
  punto 1.3). `lanciatore/dinamica.py` **non va modificato in questo step**:
  nessuna nuova firma, nessun parametro opzionale aggiunto alla sua
  `integra_ascesa_verticale`, per non introdurre rischio di regressione sui
  25 test dello Step 2 gia' verdi. Se in Fase A serve un evento diverso da
  quelli gia' cablati in `integra_ascesa_verticale`, `guida.py` chiama
  `solve_ivp` direttamente con `dinamica.derivate_stato` come `fun`,
  bypassando quella funzione di alto livello (che resta utile e testata per
  il suo caso d'uso originale, l'ascesa verticale pura completa fino a fine
  propellente).
- **implementa da zero, solo qui**, le equazioni 2D del punto 1.2
  (`derivate_stato_gravity_turn`), perche' sono fisicamente nuove (moto 2D,
  spinta proiettata lungo la velocita' via `gamma`) e non una variante
  parametrica di cio' che gia' esiste in `dinamica.py`;
- **implementa da zero, solo qui**, gli eventi terminali per lo stato a 5
  componenti della Fase B (`evento_fine_propellente_2d`:
  `y[4] - m_vuoto`, `direction=-1`; `evento_impatto_suolo_2d`: `y[1]`,
  `direction=-1`; entrambi `terminal=True`) — necessariamente reimplementati
  (non importabili da `dinamica.py`) perche' l'indice della massa/quota nel
  vettore di stato e' diverso (5 componenti vs 3), ma **concettualmente
  identici** al pattern gia' stabilito nello Step 2: questa non e'
  duplicazione della logica di drag/gravita'/massa (che resta unica, in
  `atmosfera.py`/`gravita.py`/dentro le nuove derivate), e' la
  replicazione minima e necessaria di un pattern di 1-2 righe per un
  vettore di stato diverso;
  raccomandato anche un terzo evento difensivo,
  `evento_traiettoria_invalida` (`gamma`, `direction=-1`, `terminal=True`,
  cioe' si attiva se `gamma` scende a 0 o sotto): utile proprio perche' il
  punto 1.3 dimostra che l'equilibrio a 90° e' instabile — un kick angle
  troppo grande o parametri veicolo poco realistici potrebbero far
  collassare `gamma` verso il basso (rotazione eccessiva, traiettoria che
  punta sempre piu' verso l'orizzonte fino a "tuffarsi"): questo evento
  intercetta il caso e lo rende un fallimento esplicito e diagnosticabile
  invece di un'integrazione silenziosamente patologica.

**Nota implementativa per il coder (per evitare un errore concreto):**
`solve_ivp` passa la stessa tupla `args` sia alla funzione `fun` sia a
tutte le funzioni `events`. La funzione evento "raggiunta velocita' di
kick" (Fase A) ha bisogno del parametro extra `v_kick`, che pero' non fa
parte della tupla `args` condivisa con `derivate_stato` (Step 2, firma
gia' fissata). Soluzione: legare `v_kick` con una chiusura (funzione
locale) o `functools.partial` **prima** di passare la funzione evento a
`solve_ivp`, mantenendo `terminal`/`direction` impostati come attributi
sulla funzione effettivamente passata (non sulla funzione "template").

#### 2. Struttura dei file da creare

```
lanciatore/
└── guida.py               # nuovo modulo:
                            #   - derivate_stato_gravity_turn(t, y, spinta, isp,
                            #     cd, area, mdot) — equazioni del punto 1.2
                            #   - evento_fine_propellente_2d, evento_impatto_suolo_2d,
                            #     evento_traiettoria_invalida (Fase B)
                            #   - evento_velocita_kick(...) locale per la Fase A
                            #     (chiusura/partial su v_kick, vedi nota 1.6)
                            #   - integra_gravity_turn(m0, m_vuoto, spinta, isp, cd,
                            #     area, v_kick, kick_angle_deg, h0=0.0, t_max=...)
                            #     orchestratore: Fase A (riusa derivate_stato di
                            #     dinamica.py) -> kick -> Fase B (nuove equazioni),
                            #     ritorna entrambi i risultati grezzi di solve_ivp
                            #     (nessun risultato nascosto, stesso principio
                            #     dello Step 2)

tests/
└── test_guida.py           # criteri di verifica, vedi punto 3 sotto
```

Nessuna modifica a `costanti.py`, `atmosfera.py`, `gravita.py`,
`dinamica.py` in questo step (tutti riusati cosi' come sono, vedi punto
1.6) — da verificare esplicitamente col `git diff`/review del
critic-ingegnere: zero righe modificate in quei quattro file.

#### 3. Criteri di verifica (`tests/test_guida.py`)

1. **Instabilita' dell'equilibrio a gamma=90°:** valutare
   `derivate_stato_gravity_turn` con `gamma=90°` (in radianti) e `v>0`
   qualsiasi: la componente `dgamma/dt` del vettore derivate deve essere
   `~0` (entro tolleranza numerica di macchina) — verifica diretta e
   analitica (da `cos(90°)=0`) della motivazione del kick angle (punto 1.3),
   indipendente da qualunque integrazione.
2. **Check primario (punto 1.4):** con `spinta=0`, `cd=0` (o `area=0`),
   condizioni iniziali arbitrarie sensate, integrare la Fase B per un tratto
   e verificare che `v*cos(gamma)` resti costante entro `rtol`/`atol` di
   `solve_ivp` (non una tolleranza fisica allargata).
3. **Check secondario (punto 1.4):** stesso scenario del punto 2, verificare
   che `E = v^2/2 - MU_TERRA/(R_TERRA+h)` resti costante entro la stessa
   tolleranza.
4. **Applicazione del kick:** verificare che, subito dopo il kick, lo stato
   di partenza della Fase B abbia `gamma0 == radians(90 - kick_angle_deg)`
   esatto (e' un'assegnazione diretta, non il risultato di
   un'integrazione — test esatto, non approssimato), `x0 == 0.0` esatto, e
   `h0`/`m0` **continui** con i valori finali della Fase A (uguaglianza
   esatta, nessun salto per costruzione tranne su `gamma`).
5. **Pattern difensivo `evento_traiettoria_invalida`:** con un
   `kick_angle_deg` deliberatamente esagerato (es. 60°, molto oltre il
   range realistico 1-3° del punto 1.3), verificare che l'evento si attivi
   (traiettoria che collassa verso `gamma<=0`) — analogo al test
   dell'evento "impatto suolo" dello Step 2 (verifica che il pattern
   difensivo funzioni quando serve, non solo che resti silente nel caso
   nominale).
6. **Caso nominale con parametri di test (v_kick=50, kick_angle=2°,
   parametri veicolo riusati dallo Step 2 — m0=50000 kg, m_vuoto=5000 kg,
   spinta=800000 N, isp=300 s, cd=CD da costanti.py, area=7.07 m^2):**
   nessun NaN/inf su nessuna componente di stato di entrambe le fasi;
   `m(t)` continua e coerente tra le due fasi; l'evento che ferma la Fase B
   deve essere "fine propellente" (non uno dei due eventi difensivi) — se
   uno dei due eventi difensivi scattasse con questi parametri "ragionevoli"
   e' un risultato da investigare (bug o combinazione di parametri non
   fisica), non da nascondere, coerente con la regola generale di CLAUDE.md.
7. **Riuso verificato (nota per il critic-ingegnere, non un test pytest):**
   confermare per lettura del codice che `derivate_stato_gravity_turn`
   importa `densita` e `accelerazione_gravita` (non le reimplementa), e che
   la Fase A importa/chiama `dinamica.derivate_stato` direttamente (non una
   copia).

**Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo step**
(riservato allo Step 6, come per lo Step 2) — non e' ancora presente ne'
multistadio ne' guida esoatmosferica ne' un caso di validazione con dati
reali.

#### 4. Task granulari per il coder (ordine consigliato)

1. Creare `lanciatore/guida.py`: prima `derivate_stato_gravity_turn`
   (equazioni punto 1.2, con docstring che cita Culler & Fried 1957 per
   esteso, coerentemente con lo stile gia' usato per Curtis/Anderson negli
   step precedenti).
2. Aggiungere i tre eventi della Fase B (`evento_fine_propellente_2d`,
   `evento_impatto_suolo_2d`, `evento_traiettoria_invalida`) con
   `terminal`/`direction` impostati come attributi, stesso stile dello
   Step 2.
3. Aggiungere l'evento "raggiunta velocita' di kick" per la Fase A (chiusura
   o `functools.partial`, vedi nota implementativa punto 1.6).
4. Implementare `integra_gravity_turn(...)`: Fase A con
   `dinamica.derivate_stato` riusato + nuovo evento; applicazione del kick;
   Fase B con `derivate_stato_gravity_turn` + i tre eventi; ritorno di
   entrambi i risultati grezzi di `solve_ivp` (struttura a scelta del coder,
   es. dict con chiavi `"fase_verticale"` e `"gravity_turn"`, purche'
   entrambi accessibili senza post-processing nascosto).
5. Creare `tests/test_guida.py` con i 6 criteri numerici del punto 3
   (il punto 7 e' per il critic-ingegnere, non un test pytest), parametri
   di guida marcati `# PROVVISORIO` come da punto 1.3.
6. Eseguire `pytest`, riportare solo il riepilogo pass/fail (non il log
   numerico verboso), come da regola di lavoro CLAUDE.md. Confermare che i
   25 test degli Step 1-2 restino verdi (nessuna regressione, coerente col
   punto 1.6: zero modifiche ai file esistenti).

#### 5. Note per il critic-ingegnere

- Verificare per iscritto nel codice (docstring di `guida.py`, non solo qui)
  la citazione di Culler & Fried 1957 per l'equazione di guida — vincolo
  esplicito di CLAUDE.md ("Guida, fase atmosferica ... riferimento Culler
  et al. 1957").
- Verificare che `dinamica.py`, `gravita.py`, `atmosfera.py`, `costanti.py`
  risultino **non modificati** (diff vuoto) — e' un requisito esplicito di
  questo piano (punto 1.6), non solo un auspicio.
- Verificare che il check di conservazione (punto 1.4/3.2/3.3) usi
  effettivamente `rtol`/`atol` dell'integratore come criterio di tolleranza
  e non una tolleranza fisica arbitraria "allargata per farlo passare" —
  coerente con la regola CLAUDE.md sui risultati di verifica numerica.
- Verificare che `v_kick`/`kick_angle_deg` NON siano finiti in
  `costanti.py` (sono parametri di design della guida, non costanti
  fisiche di progetto — vedi motivazione punto 1.3).
- Vincolo CLAUDE.md pertinente a questo step: "Guida, fase atmosferica:
  gravity turn (allineamento spinta-velocita' dopo un kick angle iniziale)"
  — verificare che l'angolo di attacco sia effettivamente nullo per
  costruzione nelle equazioni (nessuna componente di portanza/forza
  normale introdotta).
- **Nessun check di delta-v vs benchmark 9.1-10.0 km/s e' pertinente in
  questo step** (riservato allo Step 6) — non segnalarlo come mancante.
- Punto di attenzione esplicito: se il test nominale (criterio 6) mostra
  che l'evento `evento_traiettoria_invalida` scatta con i parametri
  "ragionevoli" proposti (kick_angle=2°), e' un segnale che va investigato
  (fisica dei parametri veicolo Step 2 non adatta al gravity turn? kick
  troppo aggressivo?) prima di rilassare qualsiasi soglia — stesso
  principio gia' seguito per il test di monotonia dello Step 2.

#### 6. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **approvabile con modifiche minori**. Segni delle
equazioni e derivazione algebrica dei due check di conservazione
verificati indipendentemente dal reviewer (corretti). Modifiche recepite:

1. **Fase A non deve bypassare silenziosamente le protezioni dello Step 2.**
   La chiamata diretta a `solve_ivp` con `dinamica.derivate_stato` in Fase A
   deve: (a) eseguire lo stesso controllo di liftoff `spinta <= m0*g(0)` →
   `ValueError` gia' presente in `integra_ascesa_verticale` (riusarlo, non
   duplicarlo silenziosamente — es. richiamare la validazione o factorizzare
   il controllo); (b) includere ANCHE `dinamica.evento_fine_propellente` e
   `dinamica.evento_impatto_suolo` tra gli eventi della Fase A, insieme al
   nuovo `evento_velocita_kick` (costo zero, gia' importabili); (c)
   `integra_gravity_turn` deve sollevare un errore esplicito se l'evento di
   kick non si attiva prima di `t_max` o prima degli altri due eventi (kick
   mai raggiunto = configurazione non valida per questo step, da segnalare
   esplicitamente, non da lasciar scorrere silenziosamente fino a t_max).
2. **Singolarita' 1/v in dgamma/dt — safeguard esplicito richiesto.** A
   differenza dello Step 2, `dgamma/dt` contiene un termine `1/v` non
   protetto. Aggiungere un evento difensivo su velocita' minima (es.
   `evento_velocita_minima`, si attiva se `v` scende sotto una soglia
   piccola ma non nulla, `terminal=True`), motivato nel codice come
   protezione dalla singolarita' analoga a `evento_traiettoria_invalida` per
   `gamma`. Estendere anche il criterio 5 (test con kick_angle=60°) con un
   controllo esplicito `np.isfinite` su tutta la traiettoria prodotta (non
   solo nel criterio 6 nominale), cosi' un collasso numerico prima
   dell'attivazione dell'evento non passa inosservato.
3. **Nuovo criterio 7 (test pytest, non solo nota per il critic-ingegnere):**
   verifica algebrica diretta che valutando `derivate_stato_gravity_turn`
   con `spinta != 0` e `drag != 0` (non il caso degenere T=0/D=0) a uno
   stato arbitrario, la componente restituita `dgamma/dt` coincida
   ESATTAMENTE con `-(g(h)/v)*cos(gamma)`, indipendente dai valori di
   spinta/drag passati — intercetta un eventuale bug che introducesse
   impropriamente un termine di spinta/drag in dgamma/dt (violerebbe il
   vincolo CLAUDE.md "angolo di attacco nullo" senza che i check di
   conservazione nel caso degenere se ne accorgano).
4. **Convenzione v0 al kick da rendere non ambigua nel codice.** Al momento
   del kick, `v0` per la Fase B deve essere preso dal valore EFFETTIVO dello
   stato all'evento (`y[1]` di Fase A al momento dell'evento), non dal
   valore letterale del parametro `v_kick` — stessa fonte (stato
   dell'evento) gia' usata per `h0`/`m0`, per coerenza ed evitare un
   disallineamento tra la tolleranza del root-finder di `solve_ivp` e il
   parametro nominale. Commentarlo esplicitamente nel codice.

Nota aggiuntiva del reviewer (non richiede modifica, solo documentazione):
il check secondario di conservazione dell'energia tratta implicitamente `h`
come distanza radiale dal centro Terra ai fini del campo gravitazionale
(gravita' sempre "verticale" nel senso locale), mentre `x` resta coordinata
piatta non curva — coerente con l'ipotesi 2D non-rotante di CLAUDE.md, ma e'
un'assunzione aggiuntiva rispetto allo Step 2 (dove il moto puramente
radiale la rendeva irrilevante). Da menzionare esplicitamente nella
docstring di `guida.py`, non solo qui.

Nota aggiuntiva del reviewer (nessuna modifica richiesta, verificato senza
problemi): nessun conflitto tra `evento_traiettoria_invalida` (gamma<=0) e
`evento_impatto_suolo` (h<=0) — dato che dh/dt=v*sin(gamma), finche' gamma
resta in (0°,90°) la quota non puo' scendere, quindi gamma deve
necessariamente attraversare lo zero prima che h possa iniziare a
diminuire. Da rendere esplicito nella docstring per un lettore futuro.

**Prossimo:** passare al sub-agente coder per l'implementazione secondo
questo piano + addendum. Al termine, aggiornare STATUS.md con l'esito
(coder + critic-ingegnere + optimizer), come per gli step precedenti.

#### 7. Esito ciclo (coder + critic-ingegnere + optimizer)

- **coder**: implementato `lanciatore/guida.py` (derivate_stato_gravity_turn,
  4 eventi Fase B, evento_velocita_kick Fase A, integra_gravity_turn) e
  `tests/test_guida.py` (7 criteri). Nessuna modifica a costanti.py/
  atmosfera.py/gravita.py/dinamica.py. 32/32 test verdi (13+12+7). Tre
  ambiguita' del piano risolte esplicitamente e segnalate (functools.partial
  per m_vuoto in evento_fine_propellente_2d; ripetizione consapevole del
  controllo di liftoff invece di importarlo, perche' dinamica.py non va
  modificato; v_min=1.0 m/s come default non specificato nel piano).
- **critic-ingegnere**: pytest rieseguito in modo indipendente (32 passed).
  Vincolo "Guida, fase atmosferica: gravity turn, riferimento Culler et al.
  1957" **verificato** — fonte citata nel codice, angolo di attacco nullo
  per costruzione, confermato anche dal test 7 (dgamma/dt indipendente da
  spinta/drag). Protezioni Fase A (liftoff, eventi difensivi Step 2)
  confermate reintrodotte, non bypassate. Safeguard 1/v verificato. Le tre
  ambiguita' risolte dal coder giudicate accettabili e ben motivate.
  Nessuna violazione. Nota: repo senza commit, non verificabile via
  git diff che i 4 moduli Step 1-2 siano invariati — verificato per lettura
  diretta (esito positivo); raccomandato un primo commit per rendere questo
  controllo automatico dai prossimi step (deciso dall'utente, non fatto
  autonomamente).
- **optimizer**: una ridondanza corretta in test_guida.py (setup duplicato
  tra test 2 e 3, estratto in costanti module-level). 32 test confermati
  verdi.

**Step 3: COMPLETATO.** Prossimo step proposto: Step 4 (guida a tangente
lineare, fase esoatmosferica).

---

### 2026-08-16 — Ciclo 4 (piano scritto direttamente dall'orchestratore, non dal planner)

**Nota di processo:** questo ciclo salta l'invocazione del sub-agente
`planner` (e, a fine ciclo, anche `optimizer`/`reporter`) — pipeline
snellita per evitare gli stalli di sessione gia' incontrati nel progetto
del collasso stellare su uno step altrettanto delicato (vedi memoria
`feedback_streamlined_pipeline.md`). Il rigore richiesto NON e' ridotto:
`reviewer` e `critic-ingegnere` restano entrambi invocati, ciascuno con
l'obbligo di ricalcolo numerico esplicito dei coefficienti A/B (vedi sotto),
non un giudizio a occhio.

**Step:** 4 — Guida a tangente lineare (fase esoatmosferica).

#### 1. Derivazione (propria, nessun accesso web disponibile ai sub-agenti)

Problema di controllo ottimo: minimizzare il tempo per raggiungere una
velocita' terminale (vx_f, vh_f) assegnata, con modulo di spinta costante,
nessun drag (fase esoatmosferica), gravita' **linearizzata a valore
costante** sull'arco di manovra (approssimazione intrinseca alla tecnica
stessa, non una scorciatoia aggiuntiva — e' esattamente cio' che rende il
risultato "lineare" invece che accoppiato a g(h) variabile).

Stato per il principio del massimo di Pontryagin: posizione (x,h) e
velocita' (vx,vh). Hamiltoniano H = 1 + λx·vx + λh·vh + λvx·a(t)cosθ +
λvh·(a(t)sinθ - g), con a(t) = spinta/m(t).

Equazioni dei costati: dλx/dt = dλh/dt = 0 (λx, λh costanti);
dλvx/dt = -λx (costante) → λvx(t) lineare in t; dλvh/dt = -λh (costante) →
λvh(t) lineare in t. Il controllo ottimo massimizza λvx·cosθ + λvh·sinθ,
quindi (cosθ,sinθ) parallelo a (λvx(t), λvh(t)): tanθ(t) = λvh(t)/λvx(t).

Scegliendo il sistema di riferimento in modo che λx=0 (asse orizzontale
lungo cui la posizione finale e' libera — scelta standard nella
derivazione classica, qui e nella letteratura PEG/Perkins), λvx(t) resta
COSTANTE e λvh(t) resta LINEARE in t, dando esattamente:

```
tan(theta(t)) = A + B*t
```

con A = λvh(0)/λvx e B = -λh/λvx, due costanti libere determinate dalle
condizioni al contorno sulla velocita' finale. **Fonte della tecnica:**
Perkins, F.M., "Derivation of Linear-Tangent Steering Laws" (citata in
CLAUDE.md) — derivazione qui rifatta da principi primi (calcolo delle
variazioni/Pontryagin), non riprodotta letteralmente dal testo originale
(nessun accesso web disponibile ai sub-agenti, stesso limite gia' gestito
per Culler & Fried allo Step 3).

**Convenzione dell'angolo:** θ e' l'angolo della SPINTA rispetto
all'orizzontale locale, in un riferimento CARTESIANO fisso (non relativo
alla velocita'). Questo e' concettualmente diverso da γ dello Step 3
(gravity turn), che era l'angolo della VELOCITA' rispetto all'orizzonte,
con spinta vincolata ad essere allineata alla velocita' stessa (angolo di
attacco nullo). Qui, al contrario, la spinta puo' NON essere allineata
alla velocita' — e' guida attiva vera, non piu' un vincolo di assetto.

#### 2. Confine di questo step (scope deciso esplicitamente)

Questo step implementa e verifica la LEGGE di guida (la dinamica sotto
tan(theta)=A+B*t con A,B DATI/noti), non la soluzione del problema al
contorno generale (trovare A,B per centrare un target orbitale reale — in
letteratura PEG questo e' tipicamente un problema accoppiato risolto per
iterazione/shooting, dato che con deplezione di massa reale l'accelerazione
a(t) non e' costante e l'integrale non ha piu' forma chiusa). Trovare A,B
per un target di missione reale resta esplicitamente FUORI SCOPE qui
(rimandato a quando servira' davvero, Step 5/6) — stessa disciplina di
scope gia' applicata allo Step 3 (nessuna anticipazione della logica di
transizione). Questo riduce il rischio di dover fare debug di un
root-finder in questo ciclo, coerente con l'obiettivo di ridurre round di
verifica per uno step gia' delicato di suo.

#### 3. Equazioni implementate

Stato `y = [x, h, vx, vh, m]` (cartesiano, non polare — qui non serve
ricostruire l'angolo della velocita', il controllo θ(t) e' gia' una
funzione esplicita ed elementare di t):

```
dx/dt  = vx
dh/dt  = vh
dvx/dt = (spinta/m) * cos(theta(t))
dvh/dt = (spinta/m) * sin(theta(t)) - g_costante
dm/dt  = -mdot
```

con `theta(t) = arctan(A + B*t)`, `g_costante` un parametro passato
esplicitamente (valutato una volta, es. con
`gravita.accelerazione_gravita(h_iniziale_fase)`, riusata per coerenza ma
NON richiamata dentro l'integrazione ad ogni passo — la costanza di g e'
l'approssimazione dichiarata di questa tecnica, vedi punto 1). Nessun
riuso diretto di `guida.derivate_stato_gravity_turn` (fisica diversa: la'
spinta e velocita' sono vincolate allineate, qui no) — nuova funzione,
stesso principio gia' seguito allo Step 3 per non riusare equazioni
fisicamente diverse solo per "risparmiare" una funzione.

#### 4. Caso di test numerico per il confronto A/B (OBBLIGATORIO,
richiesto esplicitamente dall'utente)

Con `mdot = 0` (nessuna deplezione di massa: caso limite scelto apposta
perche' rende `a(t) = spinta/m0` ESATTAMENTE costante, quindi l'integrale
ha forma chiusa esatta — non un'approssimazione), gli integrali di
`cos(theta(t))` e `sin(theta(t))` con `theta(t)=arctan(A+Bt)` si risolvono
esattamente:

```
vx(tf) - vx0 = (a0/B) * [asinh(A+B*tf) - asinh(A)]
vh(tf) - vh0 = (a0/B) * [sqrt(1+(A+B*tf)^2) - sqrt(1+A^2)] - g*tf
```

(con `a0 = spinta/m0` costante). Questa e' matematica pura (integrale di
1/sqrt(1+u^2) e u/sqrt(1+u^2)), non richiede fonte esterna — verificabile
da chiunque con carta e penna/calcolatrice, stesso principio dei check di
conservazione dello Step 3.

**Valori di test proposti** (marcare `# PROVVISORIO` in
tests/test_guida_esoatmosferica.py, come i parametri veicolo degli step
precedenti):
- `A = 0.3`, `B = -0.002 (1/s)`, `tf = 20 s`
- `m0 = 1000 kg`, `spinta = 40000 N` → `a0 = 40 m/s^2`
- `g_costante = G0 = 9.80665 m/s^2` (riuso della costante di progetto)
- `vx0 = vh0 = 0`, `x0 = h0 = 0`

**Valori CORRETTI (vedi addendum reviewer punto 8 sotto — il calcolo a mano
originale dell'orchestratore aveva un errore di arrotondamento amplificato
da cancellazione catastrofica, corretto qui, confermato indipendentemente
sia dal reviewer sia da me con `numpy.arcsinh`/`numpy.sqrt` a piena
precisione):**
- `asinh(0.3) = 0.2956730475634224`, `asinh(0.26) = 0.2571563485130149`
  (dove `0.26 = A+B*tf`)
  → `vx(tf) = (40/-0.002)*(0.2571563485130149-0.2956730475634224)
  = 770.33398100815 m/s` (NON 770.5 come da primo calcolo a mano, errore
  dovuto a cancellazione catastrofica nella sottrazione di due asinh
  vicini moltiplicata per il fattore grande a0/B=-20000 — vedi addendum)
- `sqrt(1+0.3^2) = 1.044030650891055`, `sqrt(1+0.26^2) = 1.0332473082471592`
  → `vh(tf) = (40/-0.002)*(1.0332473082471592-1.044030650891055) -
  9.80665*20 = 19.533852877918207 m/s`

**IMPORTANTE per il coder:** il test pytest (criterio 1) deve calcolare la
formula chiusa con `numpy.arcsinh`/`numpy.sqrt` IN CODICE al momento del
test, MAI hardcodare questi valori decimali arrotondati come target di
`pytest.approx` — la cancellazione catastrofica appena dimostrata prova che
valori scritti a mano con poche cifre non sono affidabili come "verita'"
hardcoded in un assert.

**Tolleranza richiesta per il confronto:** il valore calcolato dalla
formula chiusa deve coincidere con l'output di `solve_ivp` (integrazione
numerica della stessa dinamica con `mdot=0`) entro la tolleranza
dell'integratore stesso (es. `rtol=1e-8` usato per l'integrazione,
confronto con `pytest.approx(..., rel=1e-6)` per margine di accumulo,
stesso principio gia' giudicato adeguato dal critic-ingegnere allo
Step 3) — NON una tolleranza fisica allargata.

**Verifica di plausibilita' fisica del caso di test (fatta da me, per
trasparenza):** con questi parametri `dvh/dt` resta positivo per tutta la
durata (parte da `40*sin(16.7°)-9.80665 ≈ 1.68 > 0` e si riduce
gradualmente restando positivo fino a `t=tf`), quindi `h(t)` e' monotona
crescente per l'intero test — nessuna sorpresa di segno da investigare.

#### 5. Struttura dei file

```
lanciatore/
└── guida_esoatmosferica.py   # derivate_stato_tangente_lineare(t, y, A, B,
                               #   spinta, mdot, g_costante) — equazioni punto 3
                               # integra_tangente_lineare(m0, spinta, mdot, A, B,
                               #   g_costante, tf, x0=0, h0=0, vx0=0, vh0=0)
                               #   — solve_ivp su [0, tf], nessun evento
                               #   terminale (tf fisso, nessuno staging/
                               #   terminazione anticipata in questo step,
                               #   coerente col confine del punto 2)

tests/
└── test_guida_esoatmosferica.py   # criteri al punto 6
```

Nessuna modifica a `costanti.py`/`atmosfera.py`/`gravita.py`/`dinamica.py`/
`guida.py` (tutti riusati dove pertinente, non alterati).

#### 6. Criteri di verifica (`tests/test_guida_esoatmosferica.py`)

1. **Confronto A/B in forma chiusa (il check centrale, vedi punto 4):**
   con `mdot=0` e i parametri di test, `vx(tf)` e `vh(tf)` calcolati da
   `integra_tangente_lineare` devono coincidere con le formule chiuse
   entro la tolleranza dichiarata.
2. **Identita' algebrica diretta:** valutando
   `derivate_stato_tangente_lineare` a un istante/stato arbitrario,
   `dvx/dt` e `dvh/dt` restituiti devono coincidere ESATTAMENTE con
   `(spinta/m)*cos(arctan(A+Bt))` e `(spinta/m)*sin(arctan(A+Bt))-g`
   (stesso stile del criterio 7 dello Step 3, verifica diretta della
   formula indipendente dall'integrazione).
3. **Monotonia di h(t) nel caso di test nominale:** verificata
   puntualmente sulla soluzione (conseguenza attesa, vedi punto 4 —
   se fallisse sarebbe un errore nei parametri di test o nell'
   implementazione, da investigare, non da nascondere).
4. **Nessun NaN/inf** su tutta la traiettoria.
5. **Massa costante quando mdot=0:** `m(t) == m0` per tutta la traiettoria
   nel caso di test (verifica banale ma diretta che `mdot=0` sia rispettato
   e non ci sia un bug che fa decrescere la massa comunque).
6. **Caso con mdot≠0 (sanity, non hand-check):** con parametri di
   deplezione realistici (es. isp=300s coerente con G0, spinta/isp→mdot
   come negli step precedenti), nessun NaN/inf, massa decrescente
   correttamente — solo un controllo di non-patologia, il confronto
   analitico esatto resta il caso `mdot=0` del criterio 1 (la forma chiusa
   non vale quando la massa depleta davvero, come dichiarato al punto 2).
   **Protezione obbligatoria (addendum reviewer punto 8):** verificare
   ESPLICITAMENTE a monte, con un commento nel test, che `mdot*tf` resti
   ben al di sotto di `m0` per i parametri scelti (margine ampio, es.
   `mdot*tf < 0.5*m0`), cosi' che `spinta/m` non rischi mai una
   singolarita' durante l'integrazione a `tf` fisso — coerente col pattern
   gia' stabilito negli Step 2/3 (eventi/controlli difensivi su massa),
   qui realizzato come vincolo verificato sui parametri di test invece di
   un evento terminale (accettabile perche' `tf` e' fisso e breve in
   questo scenario di sanity check, non una traiettoria completa fino a
   fine propellente).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo step
(riservato allo Step 6). Nessuna soluzione del problema A/B-per-target
(riservato a quando servira', vedi punto 2).

#### 7. Nota per reviewer e critic-ingegnere (requisito esplicito
dell'utente, non facoltativo)

Ricalcolare autonomamente, con carta/calcolatrice, `vx(tf)` e `vh(tf)`
dalle formule chiuse del punto 4 usando gli stessi valori numerici, e
riportare per iscritto in STATUS.md i propri numeri calcolati e il
confronto con l'output effettivo del codice (non solo "ho controllato la
formula e sembra giusta"). Il critic-ingegnere lo fa indipendentemente dal
reviewer (sul codice implementato, non sul piano) — se i due calcoli
indipendenti (reviewer sul piano, critic-ingegnere sul codice) non
coincidono tra loro entro la tolleranza dichiarata, e' un segnale da
investigare subito, non da far passare silenziosamente.

#### 8. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **approvabile con modifiche minori**. Il reviewer ha
eseguito il ricalcolo numerico obbligatorio (a mano, con serie di Taylor +
logaritmo come doppio controllo incrociato) e ha trovato un errore reale
nel calcolo a mano dell'orchestratore: `vx(tf)` era riportato come 770.5
m/s ma il valore corretto e' 770.33 m/s (scarto relativo ~2.2e-4, sopra la
tolleranza 1e-4 richiesta), causato da cancellazione catastrofica nella
sottrazione di due `asinh` vicini moltiplicata per il fattore grande
`a0/B=-20000`. **Confermato indipendentemente anche da me (orchestratore)
con `numpy.arcsinh`/`numpy.sqrt` a piena precisione** (vedi valori corretti
gia' sostituiti al punto 4 sopra: vx(tf)=770.33398 m/s,
vh(tf)=19.53385 m/s). Modifiche recepite:

1. **Valori numerici corretti** al punto 4 (fatto, vedi sopra) — e nota
   esplicita per il coder che il test deve calcolare la formula chiusa in
   codice con numpy, mai hardcodare cifre arrotondate a mano come target.
2. **Protezione mancante per il criterio 6 (mdot≠0)** — colmata al punto 6
   sopra (verifica esplicita `mdot*tf << m0` nei parametri di test).
3. **Completare la derivazione teorica nella docstring del codice** (non
   solo qui in STATUS.md): (a) chiarire che `λx=0` segue dalla condizione
   di trasversalita' per posizione finale `x_f` LIBERA (non vincolata),
   non da una generica "scelta di sistema di riferimento" — e' una
   conseguenza del problema al contorno, non una scelta arbitraria di
   assi; (b) rendere esplicito che, per lo stesso ragionamento, se anche
   `h_f` fosse libera si avrebbe `λh≡0` e quindi `B=0`, collassando la
   legge a un angolo costante — la tangente lineare vera (B≠0) richiede
   quindi che la quota finale `h_f` sia VINCOLATA (target di iniezione
   orbitale specificato), mentre la distanza a terra `x_f` resta libera.
   Questa asimmetria (x_f libera, h_f fissata) e' l'assunzione che rende
   la derivazione autoconsistente e va scritta esplicitamente nella
   docstring di `guida_esoatmosferica.py`, non lasciata implicita.

Punto verificato senza problemi dal reviewer: convenzione dell'angolo
theta chiara; nota aggiuntiva (non blocca nulla, da menzionare in
docstring): `arctan(A+Bt)` restituisce sempre valori in (-90°,90°), quindi
la parametrizzazione impone implicitamente spinta sempre con componente
orizzontale non negativa (mai spinta "all'indietro") — ragionevole per
iniezione orbitale, ma e' un vincolo implicito della parametrizzazione da
rendere esplicito in un commento.

**Prossimo:** passare al sub-agente coder per l'implementazione secondo
questo piano + addendum.

#### 9. Esito ciclo (coder + critic-ingegnere, pipeline snellita — vedi nota di processo in testa a questo ciclo; pulizia e chiusura fatte direttamente dall'orchestratore)

- **coder**: implementati `lanciatore/guida_esoatmosferica.py`
  (derivate_stato_tangente_lineare, integra_tangente_lineare, con
  derivazione Pontryagin/costati completa nel docstring incluso
  l'addendum del reviewer) e `tests/test_guida_esoatmosferica.py` (6
  criteri, forma chiusa calcolata in codice con numpy, mai hardcodata).
  Nessuna modifica ai 5 moduli riusati (diff vuoto confermato). 38/38 test
  verdi (32+6). Criterio 1: vx(tf)=770.3339810076213 vs atteso
  770.33398100815 (scarto ~6.9e-13); vh(tf)=19.53385287791232 vs atteso
  19.533852877918207 (scarto ~3.0e-13) — ben entro rel=1e-6. Nessuna
  ambiguita' residua.
- **critic-ingegnere**: ricalcolo numerico INDIPENDENTE (dal reviewer, sul
  codice implementato, non sul piano) di vx(tf)/vh(tf) — coincidenza
  esatta con i valori del reviewer e con l'output del codice, nessuna
  discrepanza tra i due calcoli indipendenti. Vincolo "Guida, fase
  esoatmosferica: tangente lineare, rif. Perkins" **verificato** — fonte
  citata, derivazione completa incluse entrambe le parti dell'addendum
  (trasversalita' λx=0, necessita' di h_f vincolata per B≠0). Protezione
  mdot*tf<0.5*m0 del criterio 6 verificata (271.9 kg < 500 kg, margine
  ampio). pytest rieseguito in modo indipendente (38 passed). File
  riusati confermati non modificati via git diff. Nessuna violazione.
- **Pulizia (io, non optimizer):** codice e test riletti, gia' puliti e
  ben organizzati (docstring esplicative, nessuna ridondanza tra i 6
  test) — nessuna correzione necessaria, coerente con l'esito frequente
  dell'optimizer negli step precedenti.

**Step 4: COMPLETATO.** Questo era lo step piu' delicato del progetto
(coefficienti A/B della tangente lineare): la verifica rafforzata
richiesta esplicitamente dall'utente (doppio ricalcolo numerico
indipendente, reviewer sul piano e critic-ingegnere sul codice) ha
effettivamente trovato e corretto un errore reale (cancellazione
catastrofica nel primo calcolo a mano dell'orchestratore, 770.5 invece di
770.33 m/s) prima che arrivasse al codice. Prossimo step proposto: Step 5
(multistadio, eventi di staging).

---

### 2026-08-16 — Ciclo 5 (piano scritto direttamente dall'orchestratore)

**Nota di processo:** pipeline snellita confermata (vedi memoria
`feedback_streamlined_pipeline.md`): scrivo io il piano, salto
`planner`/`optimizer`/`reporter`, mantengo `reviewer` e `critic-ingegnere`
come chiamate sub-agente vere.

**Step:** 5 — Multistadio, fase atmosferica (eventi di staging, cambio
massa discontinuo).

**Decisione dell'utente sullo scope (2026-08-16, non riaprire in questo
ciclo):** lo staging qui riguarda SOLO la fase atmosferica (verticale +
gravity turn, Step 2/3), non la tangente lineare (Step 4) — quella
estensione e' stata esplicitamente aggiunta come nuovo Step 7 in
STATUS.md (non opzionale, solo rimandata), perche' un cambio di massa a
meta' della fase B renderebbe i coefficienti A/B del segmento in corso
invalidi, e ricalcolarli e' accoppiato al problema inverso gia' rimandato
allo Step 6. Implementato pero' come meccanismo GENERICO e riusabile,
cosi' la stessa struttura serve anche allo Step 7 quando ci si arrivera'.

#### 1. Design tecnico

**Rappresentazione di uno stadio:** dict con `m_prop` (kg, propellente),
`m_strut` (kg, massa strutturale espulsa al burnout), `spinta` (N), `isp`
(s). Nessuna massa "payload" esplicita — implicita in `m0` totale meno
tutto cio' che verra' bruciato/espulso.

**Riuso massimo (nessuna reimplementazione della fisica):**
- Stadio 1: riusa **direttamente e senza modifiche**
  `guida.integra_gravity_turn` (Fase A verticale + kick + Fase B gravity
  turn, Step 3), passando come `m_vuoto` la soglia di burnout DI QUESTO
  STADIO (`m0 - m_prop_1`, non la massa a vuoto finale del razzo — la
  funzione non ha bisogno di sapere che ci sono altri stadi attaccati).
- Stadi 2..N: riusano **direttamente**
  `guida.derivate_stato_gravity_turn` e i suoi eventi difensivi
  (`evento_impatto_suolo_2d`, `evento_traiettoria_invalida`,
  `evento_velocita_minima`), con un nuovo `evento_fine_propellente_2d`
  parametrizzato via `functools.partial` (stesso pattern gia' stabilito
  allo Step 3) sulla soglia di burnout di QUEL segmento.
- **Nessuna modifica a `guida.py`** (diff vuoto da verificare, come per
  ogni step precedente).

**Discontinuita' allo staging:** al burnout dello stadio i, `x, h, v,
gamma` restano ESATTAMENTE continui (solo la massa cambia); la massa
scende istantaneamente di `m_strut_i` PRIMA di iniziare l'integrazione
dello stadio i+1. Anche dopo l'ultimo stadio l'espulsione della struttura
avviene (fisicamente reale), semplicemente non c'e' integrazione
successiva.

**Nuovo modulo:** `lanciatore/staging.py`, funzione
`integra_multistadio_gravity_turn(stadi, v_kick, kick_angle_deg, h0=0.0,
t_max=1000.0, v_min=1.0)`: orchestratore che chiama
`guida.integra_gravity_turn` per lo stadio 1, poi itera sugli successivi
con `solve_ivp` diretto (dinamica + eventi riusati da `guida.py`),
applicando la discontinuita' di massa tra un segmento e il successivo.
Ritorna la lista di tutti i risultati grezzi di `solve_ivp` (nessun
risultato nascosto) piu' un riepilogo (masse/istanti di ogni evento di
staging).

#### 2. Struttura dei file

```
lanciatore/
└── staging.py   # integra_multistadio_gravity_turn(...)

tests/
└── test_staging.py   # criteri al punto 3
```

Nessuna modifica a `costanti.py`/`atmosfera.py`/`gravita.py`/
`dinamica.py`/`guida.py`/`guida_esoatmosferica.py`.

#### 3. Criteri di verifica (`tests/test_staging.py`)

1. **Continuita' di stato allo staging:** `x,h,v,gamma` identici (entro
   tolleranza numerica dell'integratore) tra l'ultimo punto dello stadio
   i e il primo punto dello stadio i+1; la massa deve differire
   ESATTAMENTE di `m_strut_i` (assegnazione diretta, test esatto, non
   approssimato).
2. **Budget di massa rispettato:** propellente bruciato + strutture
   espulse + massa finale rimanente deve tornare esattamente a `m0` —
   controllo di contabilita' indipendente dalla fisica del moto.
3. **Confronto con Tsiolkovsky per-stadio (caso limite analitico, stesso
   principio dei check di conservazione di Step 2/3/4):** il delta-v
   ideale totale (nessuna gravita'/drag) e' la SOMMA dei delta-v ideali
   di ciascuno stadio, `Isp_i*G0*ln(m_ignizione_i/m_burnout_i)` —
   identita' puramente algebrica sulle masse/Isp, verificabile a mano,
   nessuna fonte esterna necessaria. Verificare che questa somma sia
   calcolabile e che la velocita' finale INTEGRATA (con gravita'/drag
   reali) resti strettamente minore di questo limite superiore (stesso
   pattern del criterio Tsiolkovsky di Step 2).
4. **Ogni evento di burnout e' "fine propellente", non un evento
   difensivo:** con parametri di test ragionevoli, nessuno stadio deve
   fermarsi per impatto suolo/traiettoria invalida/velocita' minima — se
   succede e' un segnale da investigare, non da ignorare.
5. **Nessun NaN/inf** su tutta la traiettoria multistadio concatenata.
6. Parametri di test (2 stadi, valori plausibili) marcati
   `# PROVVISORIO`, riusando dove sensato i valori gia' usati per lo
   stadio 1 negli Step 2/3 (m0=50000, m_prop_1=45000 con m_strut_1 da
   dedurre, spinta_1=800000N, isp_1=300s; stadio 2 con massa/spinta
   ridotte in proporzione plausibile per uno stadio superiore).

Nessun confronto con il benchmark delta-v 9.1-10.0 km/s in questo step
(riservato al nuovo Step 6, dati reali). Nessuna soluzione del problema
A/B-per-target (fuori scope, deciso allo Step 4; riguarda comunque
tangente lineare, non pertinente qui).

**Prossimo:** passare al sub-agente reviewer per il controllo critico di
questo piano, poi al coder.

#### 4. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **da rivedere** (non solo modifiche minori — mancavano
parametri necessari alla firma). Modifiche recepite:

1. **`m0` (massa totale iniziale) mancante — aggiunto.** Nuova firma:
   `integra_multistadio_gravity_turn(m0, stadi, cd, area, v_kick,
   kick_angle_deg, h0=0.0, t_max=1000.0, v_min=1.0)`. `m0` e' la massa
   totale allo stacco (tutti gli stadi + payload implicito). Lo stadio 1
   riusa `guida.integra_gravity_turn(m0, m_vuoto=m0-m_prop_1, ...)`.
2. **`cd`/`area` mancanti — aggiunti come parametri GLOBALI** della
   funzione (non per-stadio): decisione di design esplicita, coerente con
   la semplificazione "Cd costante" gia' dichiarata a livello di progetto
   in CLAUDE.md (un solo Cd/area per l'intera simulazione, non solo per
   singolo stadio) — differenze aerodinamiche tra stadi (es. sgancio
   ogiva) sono un ulteriore raffinamento esplicitamente rimandato, non
   richiesto da questo step. Riusano `CD` da `costanti.py` e l'`area` gia'
   usata negli Step 2/3 (diametro 3.0 m).
3. **Massa di base per la sottrazione di `m_strut_i`: e' lo STATO
   EFFETTIVO, non il valore nominale.** Stesso principio gia' stabilito
   nell'addendum Step 3 punto 4 (v0 preso dallo stato effettivo
   dell'evento, non dal parametro nominale): la massa da cui sottrarre
   `m_strut_i` e' `risultato_stadio_i.y[4, -1]` (il valore che
   `solve_ivp` restituisce ESATTAMENTE all'evento terminale, entro la
   tolleranza del root-finder di `solve_ivp`), non il parametro nominale
   `m_vuoto_i` calcolato a monte. `x, h, v, gamma` sono presi dalla
   stessa fonte (stato effettivo), coerenza totale nella provenienza dei
   dati di continuita'.
4. **Nuovo criterio di verifica (analogo al criterio 1 dello Step 2):**
   la massa EFFETTIVA all'evento di fine-propellente dello stadio i deve
   coincidere con la soglia nominale `m_vuoto_i` entro l'`atol`
   dell'integratore (`pytest.approx` con tolleranza coerente) — questo e'
   il controllo che intercetta un eventuale disallineamento tra valore
   nominale ed effettivo, mancante nella prima stesura del piano.
5. **Tabella numerica esplicita per il caso di test a 2 stadi** (calcolata
   e verificata algebricamente, T/W e tempi di bruciamento controllati
   per plausibilita' fisica):

| Grandezza | Stadio 1 | Stadio 2 |
|---|---|---|
| `m_prop` | 40 000 kg | 8 000 kg |
| `m_strut` | 4 000 kg | 800 kg |
| `spinta` | 800 000 N | 150 000 N |
| `isp` | 300 s | 320 s |
| `mdot` (derivata) | 271.924 kg/s | 47.799 kg/s |
| `t_burn` (derivata) | ≈147.1 s | ≈167.4 s |
| T/W a ignizione (derivata) | 1.511 | 1.530 |

Masse di riferimento (derivate, da usare nei test):
`m0 = 54000 kg` (payload implicito = 1200 kg, mai nominato come campo a
se', solo cio' che resta a fine stadio 2); `m_vuoto_1 = m0 - m_prop_1 =
14000 kg`; `m_ignizione_2 = m_vuoto_1 - m_strut_1 = 10000 kg`;
`m_vuoto_2 = m_ignizione_2 - m_prop_2 = 2000 kg`; `m_finale = m_vuoto_2 -
m_strut_2 = 1200 kg` (torna esattamente al payload, verifica di
bookkeeping gia' fatta a mano). Delta-v ideali per-stadio (Tsiolkovsky,
verificati algebricamente): `dv_ideale_1 ≈ 3971.478 m/s`, `dv_ideale_2 ≈
5050.622 m/s`, somma ≈ `9022.100 m/s` (SOLO un limite superiore teorico
per il confronto del criterio 3, NON un valore da raggiungere — nessun
legame col benchmark 9.1-10.0 km/s dello Step 6, che e' un numero diverso
per un contesto diverso, guidato e con orbita reale).
6. **Criterio 2 (budget di massa) rietichettato esplicitamente come
   sanity-check sui dati di INPUT** (vero per costruzione sui parametri
   di progetto, non sulla simulazione) — il vero test sulla simulazione
   e' il nuovo criterio 4 dell'addendum sopra (massa effettiva vs
   nominale a ogni evento).
7. **Nota per il critic-ingegnere (punto minore del reviewer):** con piu'
   stadi, `gamma` puo' avvicinarsi alla soglia difensiva
   `evento_traiettoria_invalida` (gamma<=0) in modo piu' pronunciato che
   con un solo stadio — verificare nel caso di test nominale che questo
   NON accada (coerente col criterio 4 gia' nel piano, "ogni evento di
   burnout e' fine_propellente, non un evento difensivo").

#### 5. Scoperta empirica del coder e correzione dei parametri di test
(risultato inatteso investigato secondo la regola CLAUDE.md — NON
un aggiustamento silenzioso di tolleranze)

Con la tabella numerica dell'addendum punto 5 (stadio 2: m_prop=8000,
m_strut=800), il coder ha implementato correttamente `staging.py` (diff
vuoto sui moduli riusati, logica di continuita' conforme al piano), ma
**2 test su 6 sono falliti**: lo stadio 2 si fermava per l'evento
difensivo `evento_traiettoria_invalida` (gamma collassato a zero) invece
che per `fine_propellente`, a t≈75.7s contro un tempo di bruciamento
nominale atteso di ≈167.4s. Il coder ha correttamente RIFIUTATO di
alterare le tolleranze o i parametri per far passare i test in silenzio,
e ha segnalato il problema per decisione.

**Causa fisica identificata (verificata numericamente da me,
orchestratore):** l'equazione del gravity turn,
`dgamma/dt = -(g(h)/v)*cos(gamma)`, NON dipende dalla spinta — l'angolo
di rotta decresce monotonicamente verso zero col tempo per costruzione,
a un tasso determinato solo da `g/v` e dall'angolo corrente, qualunque
sia il livello di spinta dello stadio in corso. Concatenare due
bruciamenti lunghi sotto la STESSA legge di guida non e' automaticamente
valido: lo stadio 1 (147.1s di gravity turn) porta gamma gia' a ≈12° a
fine bruciamento (v≈2899 m/s, h≈54.4 km); un secondo bruciamento di
167.4s supplementari sotto la stessa legge fa collassare gamma sotto
zero circa a meta' strada. Non e' un bug: e' esattamente il motivo per
cui un lanciatore reale smette di usare il gravity turn ben prima di
esaurire piu' stadi consecutivi sotto la stessa legge — coerente col
fatto che il progetto passa a un'altra legge di guida (tangente lineare,
Step 4) proprio per la fase successiva. Per QUESTO step (staging isolato
dalla scelta di quando cambiare legge di guida, per scope esplicito) la
tabella di test deve semplicemente restare entro la finestra in cui il
gravity turn concatenato regge, senza introdurre alcuna logica di
cambio-legge (fuori scope).

**Ricerca numerica (fatta da me) di parametri stadio 2 auto-consistenti:**
a spinta/isp fissati (150000 N, 320 s, stessi valori dell'addendum), ho
cercato il valore di `m_prop_2` tale che il tempo di bruciamento nominale
resti sotto la soglia di collasso di gamma, con margine di sicurezza
(6500 kg regge fino a t=136.0s completando il nominale; 7000 kg fallisce
a t=140.2s contro un nominale di 146.4s — il punto di rottura e' tra i
due). **Scelto `m_prop_2 = 6000 kg`** (margine confortevole sotto 6500),
`m_strut_2 = 600 kg` (stesso rapporto 10:1 di massa strutturale/propellente
dello stadio 1).

**Tabella numerica CORRETTA (sostituisce quella del punto 4.5
dell'addendum, tutti i valori ricalcolati e verificati end-to-end con
`integra_multistadio_gravity_turn` effettivamente eseguita, non solo
algebra a parte):**

| Grandezza | Stadio 1 | Stadio 2 |
|---|---|---|
| `m_prop` | 40 000 kg | **6 000 kg** |
| `m_strut` | 4 000 kg | **600 kg** |
| `spinta` | 800 000 N | 150 000 N |
| `isp` | 300 s | 320 s |
| `mdot` (derivata) | 271.924 kg/s | 47.799 kg/s |
| `t_burn` nominale (derivata) | ≈147.1 s | ≈125.5 s |
| T/W a ignizione (derivata) | 1.575 | 1.961 |

`m0 = 51800 kg` (payload implicito = 1200 kg, invariato); `m_vuoto_1 =
m0 - m_prop_1 = 11800 kg`; `m_ignizione_2 = m_vuoto_1 - m_strut_1 = 7800
kg`; `m_vuoto_2 = m_ignizione_2 - m_prop_2 = 1800 kg`; `m_finale =
m_vuoto_2 - m_strut_2 = 1200 kg` (= payload, bookkeeping verificato).
Delta-v ideali per-stadio: `dv_ideale_1 ≈ 4352.066 m/s`, `dv_ideale_2 ≈
4601.553 m/s`, somma ≈ `8953.619 m/s` (limite superiore teorico, non un
target).

**Verifica end-to-end eseguita da me con questi valori:** entrambi gli
stadi terminano per `fine_propellente` (differenza massa
effettiva/nominale = 0.0 in entrambi i casi, entro la precisione di
macchina); stato finale: v≈7508.7 m/s, gamma≈8.10° (positivo, margine
sano rispetto alla soglia difensiva), v_finale < somma ideale (8953.6
m/s) come atteso dal criterio 3.

**Azione per il coder:** aggiornare SOLO i valori numerici in
`tests/test_staging.py` (m_prop_2, m_strut_2 e le masse derivate) secondo
questa tabella corretta — nessuna modifica a `staging.py` (la logica era
gia' corretta, il problema era solo nei parametri di test). Rieseguire
`pytest` per confermare 6/6 verdi in `test_staging.py`.

#### 6. Esito ciclo (coder x2 + critic-ingegnere, pipeline snellita —
pulizia e chiusura fatte direttamente dall'orchestratore)

- **coder (1° passaggio):** implementato `lanciatore/staging.py` e
  `tests/test_staging.py` secondo il piano + addendum. Diff vuoto sui 6
  moduli riusati confermato. 2 test su 6 falliti per un motivo fisico
  reale (gamma collassato, vedi punto 5) — il coder ha correttamente
  rifiutato di aggiustare tolleranze o parametri in silenzio e ha
  segnalato il problema, esattamente come richiesto dalla regola di
  lavoro CLAUDE.md sui risultati inattesi.
- **Investigazione causa fisica (io, orchestratore):** confermata la
  causa (dgamma/dt indipendente dalla spinta), trovata numericamente una
  tabella di parametri auto-consistente con margine (m_prop_2=6000,
  m_strut_2=600), verificata end-to-end prima di rimandare al coder.
- **coder (2° passaggio):** aggiornati solo i numeri in
  `tests/test_staging.py` secondo la tabella corretta, nessuna modifica a
  `staging.py`. 44/44 test verdi. Segnalata (senza correggerla, come da
  istruzione) un'incongruenza minore nella docstring di `staging.py`
  (sezione "Solleva RuntimeError" residua da una versione precedente non
  più valida) — corretta direttamente da me (edit di sola docstring,
  nessuna modifica di logica).
- **critic-ingegnere:** pytest rieseguito in modo indipendente (44
  passed). Vincolo "Multistadio: eventi di staging come eventi terminali
  ODE" **verificato**, implementazione pulita. Continuità di stato via
  stato EFFETTIVO (non nominale) verificata riga per riga. Tsiolkovsky
  per-stadio ricalcolato indipendentemente (coincidenza esatta). Ha
  riprodotto autonomamente sia il caso di collasso con i parametri
  originali (t=75.683s, coincide) sia il punto di rottura empirico
  (6500kg regge, 7000kg no) sia il caso corretto finale — giudicata
  genuina l'investigazione della causa fisica, non un numero scelto a
  caso. Nessuna violazione.
- **Pulizia (io, non optimizer):** codice e test riletti, già puliti,
  nessuna correzione necessaria oltre al fix di docstring sopra.

**Step 5: COMPLETATO.** Il ciclo ha prodotto un secondo esempio (dopo la
cancellazione catastrofica dello Step 4) di un risultato inatteso
genuinamente investigato invece che nascosto: qui la fisica del gravity
turn (angolo che collassa indipendentemente dalla spinta) ha reso non
validi i primi parametri di test scelti, e la correzione richiesta una
ricerca numerica esplicita del punto di rottura, non un ritocco casuale.
Prossimo step proposto: Step 6 (caso di validazione con dati reali +
confronto delta-v vs benchmark 9.1-10.0 km/s).

---

### 2026-08-16 — Ciclo 6 (piano scritto direttamente dall'orchestratore)

**Nota di processo:** pipeline snellita confermata: scrivo io il piano,
salto `planner`/`optimizer`/`reporter`, mantengo `reviewer` e
`critic-ingegnere` come chiamate sub-agente vere.

**Riorganizzazione roadmap (fatta prima di questo piano, vedi sezione
"## Step" in testa al file):** inserito nuovo Step 6 dedicato al problema
inverso della tangente lineare, scorporandolo da quello che era lo Step
6 "Caso di validazione" (ora Step 7) — coerente con la nota utente
lasciata al Ciclo 4/5 ("il planner deve trattarlo come uno step a se
stante quando ci si arriva"). Rinumerati a scendere anche Step 8
(estensione staging in tangente lineare, gia' dipendeva dal problema
inverso), 9 (visualizzazione), 10 (validazione/limiti), 11 (pulizia).

**Step (nuovo numero):** 6 — Guida a tangente lineare, problema inverso
(root-finding A/B per un target di velocita' terminale, deplezione di
massa reale).

#### 1. Design tecnico

**Problema:** dati stato iniziale (x0,h0,vx0,vh0,m0), spinta, mdot,
g_costante, tempo di volo `tf` FISSO (non un'incognita aggiuntiva,
coerente con lo scope gia' stabilito allo Step 4), e una velocita'
terminale desiderata (vx_target, vh_target), trovare A, B tali che
integrando `guida_esoatmosferica.derivate_stato_tangente_lineare`
(RIUSATA SENZA MODIFICHE) su [0,tf] si ottenga vx(tf)≈vx_target,
vh(tf)≈vh_target.

**Perche' serve root-finding:** con `mdot≠0` reale, `a(t)=spinta/m(t)`
non e' costante, quindi la forma chiusa dello Step 4 (valida solo per
`mdot=0`) non si applica — esattamente la ragione per cui lo Step 4 ha
esplicitamente rimandato questo problema (vedi la sua docstring,
"Confine di questo step").

**Metodo:** `scipy.optimize.fsolve` su residuo
`R(A,B) = (vx(tf;A,B)-vx_target, vh(tf;A,B)-vh_target)`, dove
`vx(tf;A,B), vh(tf;A,B)` si ottengono chiamando
`guida_esoatmosferica.integra_tangente_lineare` (RIUSATA SENZA
MODIFICHE) con i valori correnti di A,B. Nessuna reimplementazione della
dinamica.

**Guess iniziale:** dalla formula chiusa dello Step 4 (caso limite
mdot=0, `a0=spinta/m0`), punto di partenza fisicamente motivato invece
di 0,0 arbitrario.

**Fallimento esplicito:** se `fsolve` non converge (`ier != 1` o residuo
finale sopra una tolleranza dichiarata), sollevare `RuntimeError`
esplicito col residuo raggiunto — mai accettare silenziosamente una
soluzione non convergente.

**Nuova funzione, stesso modulo Step 4 (non un nuovo file):** aggiunta a
`lanciatore/guida_esoatmosferica.py` (le due funzioni esistenti restano
INVARIATE, solo aggiunta — stesso principio gia' usato per `CD` in
`costanti.py` allo Step 2):
`risolvi_coefficienti_tangente_lineare(m0, spinta, mdot, g_costante, tf,
vx_target, vh_target, x0=0.0, h0=0.0, vx0=0.0, vh0=0.0, A0=None,
B0=None)`.

#### 2. Verifica (progettata per essere non tautologica)

Una semplice verifica "risolvi per un target, poi integra e controlla
che il residuo sia piccolo" rischia di essere tautologica (verifica solo
il criterio di arresto del solver, non la correttezza). Test case
progettato cosi':

1. `A_vero=0.3`, `B_vero=-0.002` (stessi valori dello Step 4, per
   continuita'), parametri con `mdot≠0` REALE: `spinta=40000N,
   m0=1000kg, isp=300s` (→ `mdot≈13.61 kg/s`), `tf=20s` (margine di
   massa ampio, gia' verificato allo Step 4 criterio 6).
2. Integrazione AVANTI con A_vero,B_vero (via `integra_tangente_lineare`,
   gia' verificata) per ottenere (vx_target, vh_target) — risultato
   NUMERICO (non forma chiusa, non esiste con mdot≠0): il target e'
   garantito raggiungibile per costruzione.
3. Il solver parte da un guess iniziale DIVERSO dalla verita' (il guess
   automatico dalla formula chiusa mdot=0, vicino ma non identico ad
   A_vero/B_vero perche' ignora la deplezione di massa).
4. **Verifica non tautologica:** il solver deve RECUPERARE A,B vicini ad
   A_vero,B_vero (tolleranza esplicita, es. 1e-4 relativo) — non solo un
   residuo piccolo sul target, ma il recupero dei parametri noti che
   hanno generato il target. Verifica aggiuntiva: re-integrando con A,B
   risolti si riottiene vx_target/vh_target entro la tolleranza del
   solver.
5. Test di fallimento esplicito: un target irraggiungibile con questi
   spinta/tf (es. velocita' finale enormemente superiore a quanto la
   spinta puo' fornire in tf secondi) deve produrre `RuntimeError`, non
   un risultato silenzioso e sbagliato.

#### 3. Struttura dei file

```
lanciatore/
└── guida_esoatmosferica.py   # SOLO AGGIUNTA:
                               #   risolvi_coefficienti_tangente_lineare(...)
                               #   le 2 funzioni Step 4 restano invariate

tests/
└── test_guida_esoatmosferica.py   # nuovi test aggiunti (criteri punto 2),
                                     # i 6 test Step 4 restano invariati
```

Nessuna modifica a nessun altro modulo. Nessun confronto col benchmark
delta-v 9.1-10.0 km/s in questo step (resta al nuovo Step 7, dati
reali).

**Prossimo:** passare al sub-agente reviewer per il controllo critico di
questo piano, poi al coder.

#### 4. Addendum (reviewer + investigazione numerica indipendente
dell'orchestratore) — recepito prima di passare al coder

Verdetto reviewer: **approvabile con modifiche**. Il reviewer non aveva
accesso a un interprete Python in quella sessione, quindi ha segnalato i
suoi dubbi come rischi da verificare empiricamente, non fatti accertati
(punto 2: conditioning del problema; punto 3: ambiguita' del "guess dalla
formula chiusa"). **Li ho verificati io stesso con Bash — e ho trovato un
problema PIU' SERIO di quanto il reviewer sospettasse.**

**Scoperta (risultato inatteso, investigato secondo la regola CLAUDE.md,
non nascosto):** ho costruito il target integrando avanti con
A_vero=0.3, B_vero=-0.002 (mdot reale), poi ho fatto risolvere il
problema inverso a `scipy.optimize.fsolve` partendo da un guess
"ragionevole" (angolo medio verso il target, `B0=0`). **`fsolve` converge
con residuo a precisione di macchina, ma su una soluzione SBAGLIATA**:
`A=0.258, B=+0.002` invece di `A=0.3, B=-0.002` (segno di B INVERTITO,
errore relativo su A del 14%, su B del 200%). Ripetuto con altre coppie
A_vero/B_vero (es. 0.5/-0.02, 0.4/-0.015): lo stesso fenomeno si presenta
sistematicamente, non e' un caso isolato.

**Causa:** il problema (A,B)→(vx(tf),vh(tf)) ammette (almeno) DUE radici
distinte per gli stessi due target — una "vera" e una "spuria" con B
approssimativamente di segno opposto (una sorta di soluzione speculare:
un profilo theta(t) leggermente crescente invece che leggermente
decrescente puo' produrre pressoche' la stessa media integrata su un
arco breve). Ho verificato con un test di consistenza: **il guess
iniziale determina su quale radice converge Newton/fsolve** — un guess
con `B0` dello STESSO SEGNO della verita' (anche solo approssimato, es.
scalato 0.7x o 1.3x del valore vero) converge sempre alla radice
corretta con residuo a precisione di macchina; un guess con `B0=0` o
segno opposto converge quasi sempre alla radice spuria. **Il tentativo
di costruire un guess automatico dalla formula chiusa dello Step 4
(risolvendo il sistema implicito asinh/sqrt) NON risolve il problema**:
l'ho implementato e testato, e converge anch'esso in modo incoerente
sulla radice giusta o sbagliata a seconda del caso — non e' un problema
di qualita' del guess "quanto vicino", ma di quale BACINO di attrazione
(segno di B) si imbocca, informazione che il guess automatico non ha
modo di indovinare senza gia' conoscere la risposta.

**Conseguenza per il design (revisione del punto 1 sopra):**
- **`A0`, `B0` diventano parametri OBBLIGATORI** di
  `risolvi_coefficienti_tangente_lineare` (non piu' calcolati
  automaticamente da una formula chiusa): la funzione documenta
  esplicitamente che la convergenza al risultato fisicamente inteso
  dipende da un guess nel bacino corretto (tipicamente noto dal contesto
  di missione — es. il segno atteso della variazione dell'angolo di
  guida — o da una soluzione vicina gia' nota, come avviene nei cicli
  iterativi reali di guida esplicita: PEG non riparte mai da un default
  generico, riusa la soluzione del ciclo precedente). Questo NON e' un
  compromesso per pigrizia: e' una proprieta' strutturale del problema,
  documentata invece che nascosta.
- **Safeguard di robustezza (non risolve l'ambiguita' globale, la
  rileva):** dopo la convergenza dal guess fornito, la funzione ripete
  la risoluzione da un SECONDO guess ottenuto perturbando quello fornito
  entro lo stesso bacino (es. `A0*1.3, B0*1.3`, stesso segno) e verifica
  che le due soluzioni coincidano entro tolleranza stretta. Se non
  coincidono, `RuntimeError` esplicito (bacino instabile/ambiguo con
  questo guess, non un risultato silenzioso e potenzialmente sbagliato).
- **Tolleranza di residuo:** fissata a `1e-6` (assoluta, sulle componenti
  di velocita' in m/s) per la condizione di successo di `fsolve`
  (verificata sia su `ier==1` sia sul residuo effettivo `info['fvec']`,
  non fidandosi del solo flag, come gia' segnalato dal reviewer).
- **Test companion non circolare (punto 1 del reviewer, risolto):**
  aggiungere un test con `mdot=0` in cui il target e' calcolato con la
  formula chiusa ESATTA dello Step 4 (non con `integra_tangente_lineare`
  con mdot reale), rompendo la circolarita' tra generazione del target e
  funzione di ricerca.
- **Target irraggiungibile (punto 6 del reviewer, risolto con numero
  esplicito):** limite di Tsiolkovsky per questi parametri,
  `Δv_ideale = Isp*G0*ln(m0/m(tf)) ≈ 933 m/s` (verificato:
  `300*9.80665*ln(1000/728.076)`). Usare `vh_target` enorme (non
  `vx_target`, per evitare la regione di non-monotonicita' in A segnalata
  dal reviewer punto 6) ben oltre questo limite, es. `vh_target = 5000`
  m/s, per un fallimento di convergenza garantito e pulito.
- **Test di documentazione del fenomeno (nuovo, non nel piano
  originale):** un test dedicato che mostra ESPLICITAMENTE che un guess
  con segno di B0 sbagliato converge su una soluzione diversa da
  A_vero/B_vero pur con residuo piccolo — non per dimostrare un difetto,
  ma per verificare che il fenomeno documentato nella docstring sia
  realmente riproducibile e non un'illazione, e per proteggere contro
  una futura "correzione" involontaria che lo nasconda.

**Valori numerici verificati da me (Bash, per il coder):** con
`m0=1000, spinta=40000, isp=300, tf=20` (`mdot≈13.609 kg/s`,
`mdot*tf≈272.1 kg`, margine ampio come da Step 4), `A_vero=0.3,
B_vero=-0.002`: `vx_target=899.267274 m/s`, `vh_target=54.681127 m/s`
(calcolati da `integra_tangente_lineare`, NON a mano — nessun rischio di
cancellazione catastrofica, il valore e' quello che il codice produce
realmente). Guess raccomandato per il test: `A0=A_vero*0.7=0.21,
B0=B_vero*0.7=-0.0014` (converge esatto, verificato: residuo ~1e-13).

#### 5. Esito ciclo (coder + critic-ingegnere, pipeline snellita —
pulizia e chiusura fatte direttamente dall'orchestratore)

- **coder**: implementata `risolvi_coefficienti_tangente_lineare` in
  coda a `guida_esoatmosferica.py` (funzioni Step 4 invariate,
  confermato), con il safeguard a doppio guess e la gestione esplicita
  del fallimento esattamente come da addendum. 4 nuovi test aggiunti
  (recupero non tautologico, companion non circolare mdot=0, fallimento
  esplicito con limite Tsiolkovsky calcolato in codice, documentazione
  del fenomeno delle radici spurie). 48/48 test verdi, valori numerici
  coincidenti al bit con la mia verifica indipendente.
- **critic-ingegnere**: ha RIPRODOTTO DA ZERO (script separato, non
  copiato da STATUS.md) il fenomeno delle radici multiple — stessi
  numeri esatti (vx_target/vh_target coincidenti, radice spuria
  A=0.257901/B=+0.001999 coincidente), confermando che non è un
  artefatto di trascrizione. Verificato anche che il safeguard non è
  codice morto (testato con guess al bordo del bacino, scatta
  effettivamente). pytest rieseguito in modo indipendente (48 passed).
  Nessuna violazione. Nota: nessun commit intermedio dopo lo Step 4
  rende impossibile un `git diff` letterale per queste verifiche —
  raccomandato un commit dopo ogni step completato (proporrò ora
  all'utente di committare Step 4-6 insieme).
- **Pulizia (io, non optimizer):** codice già ben strutturato
  (funzione + helper privato separati, docstring esaustiva), nessuna
  correzione necessaria.

**Step 6: COMPLETATO.** Ciclo che ha prodotto la scoperta più
significativa del progetto finora dal punto di vista numerico: il
problema inverso della tangente lineare ammette radici multiple/spurie,
un fenomeno reale (verificato indipendentemente due volte, da me e dal
critic-ingegnere) non documentato nella letteratura di riferimento
citata in modo esplicito per questo caso — la soluzione adottata (guess
obbligatorio nel bacino corretto + safeguard a doppio guess che rileva,
senza pretendere di risolvere, l'ambiguità) è un design onesto rispetto
al problema reale, non un aggiramento. Prossimo step proposto: Step 7
(caso di validazione con dati reali, ora può riusare questo risolutore).

---

### 2026-08-16 — Ciclo 7 (piano scritto direttamente dall'orchestratore, indagine numerica pre-piano)

**Nota di processo:** pipeline snellita confermata. Prima di scrivere
questo piano ho fatto un'indagine numerica diretta (Bash, script ad-hoc)
per de-rischiare la parte più incerta dello step — vedi sotto — invece di
scoprirla solo a implementazione avviata.

**Step:** 7 — Caso di validazione con dati reali di un lanciatore
pubblico + confronto delta-v vs benchmark 9.1-10.0 km/s.

#### 0. Scoperta pre-piano: gap fisico nel modello a gravità costante

Recuperati dati pubblici (Wikipedia, "Falcon 9" Block 5, consultata
2026-08-16, https://en.wikipedia.org/wiki/Falcon_9) e tentato
l'assemblaggio Stadio1 (gravity turn) + Stadio2 (tangente lineare) PRIMA
di scrivere il piano. **Trovato:** con gravità COSTANTE (Step 4/6), la
velocità orbitale a 200 km (~7788 m/s) non è raggiungibile dallo
Stadio 2 — verificato con ricerca su griglia estesa A/B, deficit ~350
m/s anche a payload zero. **Causa:** vicino alla velocità orbitale, il
sollievo centripeto `vx²/(R+h)` diventa comparabile alla gravità stessa
(per definizione `v_orbitale²/r = g(r)`), effetto ignorato dal modello a
gravità costante.

Presentata la scoperta all'utente con `AskUserQuestion` (due opzioni:
aggiungere il termine centripeto, o tenere gravità costante e validare
solo la capacità ideale del veicolo con arithmetica pura). **L'utente ha
scelto di aggiungere il termine centripeto.**

**Verificato che la correzione risolve il problema** (non solo in
teoria): con `g_eff(h,vx) = accelerazione_gravita(h) - vx²/(R_TERRA+h)`,
e correggendo anche un bug nel mio script di verifica (mancava
l'espulsione della massa strutturale dello Stadio 1 prima dell'ignizione
dello Stadio 2 — errore di massa del 21%, non un problema del progetto),
il target si raggiunge quasi esattamente: `vx=7788.487984949658` contro
target `7788.487984973157` (scarto ~2.3e-8 m/s), `vh≈-1e-8`, a
`A≈0.2197, B≈-0.000426`, quota finale ≈195.7 km.

#### 1. Fase A — Correzione fisica (PRIMA modifica di codice già testato
nel progetto)

Fino ad ora ogni step ha SOLO aggiunto codice, mai modificato funzioni
già verificate. Prima eccezione, motivata da un gap fisico reale
verificato empiricamente, approvata dall'utente. Stesso rigore di uno
step delicato come lo Step 4 originale.

#### 1bis. Addendum (reviewer) — correzione era INCOMPLETA, recepito
prima di passare al coder

Verdetto reviewer: **da rivedere**, non modifiche minori. Il reviewer ha
derivato le equazioni polari complete e dimostrato che la mia proposta
iniziale (correggere solo `dvh/dt`) omette un termine gemello in
`dvx/dt`, violando conservazione di energia e momento angolare nel
limite T=0/D=0. **Ho verificato io stesso numericamente (Bash,
`scipy.integrate.solve_ivp` a `rtol=atol=1e-12`):** con la correzione
COMPLETA, E ed L si conservano a precisione di macchina (~1e-14/1e-15
di variazione relativa su 200s); con la correzione a meta' proposta
inizialmente, E varia dello 0.197% e L dell'1.25% sugli stessi 200s —
un errore reale, non accademico, confermato.

**Design CORRETTO e completo (sostituisce il punto 1 originale):**

```
dx/dt  = vx
dh/dt  = vh
dvx/dt = (spinta/m)*cos(theta) - vh*vx/(R_TERRA+h)     [NUOVO termine]
dvh/dt = (spinta/m)*sin(theta) - g(h) + vx**2/(R_TERRA+h)
dm/dt  = -mdot
```

con `g(h) = gravita.accelerazione_gravita(h)` (riusata, non ricalcolata).
Entrambi i termini derivano dalla stessa trasformazione in coordinate
locali (x tangenziale, h radiale) attorno a un corpo sferico — non e'
opzionale aggiungerne uno solo, sono i due lati della stessa identita'
geometrica (derivazione completa nella docstring, con la dimostrazione
di conservazione E/L come prova).

**Nuova strategia di verifica (sostituisce l'idea originale "avvicinamento
alla vecchia forma chiusa a bassa velocità", giudicata dal reviewer
insensibile proprio al tipo di bug trovato):**
- **Check di conservazione nel limite T=0/D=0** (stesso principio
  esatto di Step 3, ora applicato alle equazioni corrette): con
  `spinta=0`, verificare che `vx*(R_TERRA+h)` (momento angolare
  specifico) ed `E=(vx**2+vh**2)/2 - MU_TERRA/(R_TERRA+h)` (energia
  specifica) restino costanti entro `rtol`/`atol` dell'integratore —
  non una tolleranza fisica allargata. Sostituisce il vecchio confronto
  a forma chiusa come check primario.
- **Target per i test non circolari (Step 6) generati SEMPRE per
  integrazione avanti** (mai forma chiusa): dato che ora `dvx/dt`
  dipende anch'esso dalla quota/velocità radiale, non esiste più alcuna
  forma chiusa nemmeno nel caso `mdot=0` — unifica la strategia già
  usata per il caso `mdot≠0` di Step 6 (target noto da A_vero/B_vero
  integrati avanti, poi recupero verificato dal solver).
- La vecchia forma chiusa (asinh/sqrt) e i vecchi criteri 1/8 di
  Step 4/6 che la usavano sono SUPERATI, non "riavvicinati" — vanno
  riscritti, non adattati con una tolleranza più larga.

**Task espliciti per il coder (dettaglio richiesto dal reviewer, non
lasciato implicito):**
1. Aggiornare `derivate_stato_tangente_lineare` con ENTRAMBI i nuovi
   termini; rimuovere `g_costante` dalla firma; importare `R_TERRA` da
   `costanti` (oltre a `accelerazione_gravita` da `gravita`, già
   importato).
2. Aggiornare `integra_tangente_lineare`, `risolvi_coefficienti_tangente_lineare`
   E l'helper privato `_residuo_velocita_finale` (il reviewer ha
   verificato che anche questo chiama `g_costante` — dimenticarlo
   avrebbe rotto la catena).
3. **TUTTI e 10 i test esistenti** in `tests/test_guida_esoatmosferica.py`
   (non solo 4) passano `g_costante` esplicitamente e vanno aggiornati
   per rimuoverlo — verificato dal reviewer leggendo il file per intero.
4. `test_2_identita_algebrica_diretta`: deve ricalcolare `g_eff`
   (entrambi i termini) IN MODO INDIPENDENTE nel test, non copiare la
   formula dall'implementazione (altrimenti diventa tautologico).
5. `test_8_companion_non_circolare_...`: riprogettato secondo la nuova
   strategia (target da integrazione avanti, non forma chiusa).
6. `test_9_fallimento_esplicito_target_irraggiungibile`: resta valido
   solo di firma (il limite di Tsiolkovsky è indipendente dal modello di
   gravità per costruzione) — confermato dal reviewer, nessuna modifica
   concettuale.
7. Aggiungere i 2 nuovi test di conservazione (vx*(R+h) ed E) descritti
   sopra.
8. Verificare numericamente (non assumere per analogia) se il fenomeno
   delle radici multiple/spurie di Step 6 si ripresenta con la nuova
   fisica — la topologia del sistema residuo potrebbe essere cambiata.

**Asimmetria dichiarata col gravity turn (Step 3, NON toccato in questo
ciclo):** `guida.py` continua a usare `dgamma/dt=-(g/v)*cos(gamma)`
senza sollievo centripeto. Quantificato (non solo stimato) con i dati
reali di Step 5/Ciclo 7: a fine Stadio 1, v²/(R+h) è il **12.2%** di
g(h) con dati SL (v≈2756 m/s, h≈45.3km) e il **17.3%** con dati vuoto
(v≈3266 m/s, h≈80km stimata) — non trascurabile ma sensibilmente più
piccolo del quasi-100% raggiunto vicino alla velocità orbitale (dove il
gap era scoperto). Lasciato così per questo ciclo (nessuna modifica a
codice di Step 3 oltre lo scope approvato), dichiarato esplicitamente
come approssimazione aggiuntiva nella docstring, non nascosto.

#### 2. Fase B — Assemblaggio validazione con dati reali

**Dati veicolo (fonte citata nel codice):**

| | Stadio 1 | Stadio 2 |
|---|---|---|
| m_strut | 25 600 kg | 3 900 kg |
| m_prop | 395 700 kg | 92 670 kg |
| spinta (vuoto) | 8 227 000 N | 981 000 N |
| Isp (vuoto) | 312 s | 348 s |

Payload: 22 800 kg (capacità LEO max, stessa fonte). Diametro 3.7 m →
area ≈10.75 m² (CD=0.3 già in costanti.py).

**Perché Isp/spinta da VUOTO per entrambi gli stadi:** con valori SL per
lo Stadio 1, Δv ideale totale ≈8764 m/s, SOTTO 9.1-10.0 km/s. Con valori
da vuoto (più rappresentativi per la maggior parte del volo dello
Stadio 1, ad alta quota): `4027.40 + 5110.76 = 9138.15 m/s` — dentro il
range. `m0=540 670 kg`, `m_vuoto_1=144 970 kg`, `m_ignizione_2=119 370
kg`, `m_vuoto_2=26 700 kg` (payload finale, bookkeeping verificato).

**Target orbitale:** h=200 km, `v_orbitale = sqrt(MU_TERRA/(R_TERRA+h))
≈ 7788.49 m/s` (dalle costanti già citate, Step 1).

**Confine delle fasi di guida (deciso ora):** Stadio 1 = gravity turn
intero (riuso non modificato `guida.integra_gravity_turn`, v_kick=50,
kick_angle_deg=2.0); Stadio 2 = tangente lineare intero (equazioni
corrette Fase A), A/B dal risolutore Step 6, `tf` = bruciamento completo
Stadio 2.

**Guess iniziale (verificato):** `A0≈0.22, B0≈-0.0004`. Se non converge
per differenze minori di implementazione, esplorare a partire da
`A0=tan(gamma_fine_stadio1)`, B0 piccolo negativo — documentare quanti
tentativi servono, non nascondere.

**Nuovo modulo:** `lanciatore/validazione.py`.

#### 3. Criteri di verifica

1. Δv ideale totale (Tsiolkovsky, calcolato in codice) in 9.1-10.0 km/s
   — verificato ≈9138.15 m/s a mano.
2. Traiettoria simulata raggiunge il target entro tolleranza stretta.
3. Scomposizione perdite (numpy.trapz su g(h(t))·sin(gamma(t)) per
   gravità — gamma effettivo di entrambe le fasi; su D(t)/m(t) per drag,
   solo fase atmosferica) confrontata con 1.0-1.5 (gravità) + 0.1-0.4
   (drag) km/s attesi. Eventuale residuo attribuito esplicitamente a
   "perdita di manovra" (thrust non allineato alla velocità in tangente
   lineare) se non trascurabile — non nascosto.
4. Nessun NaN/inf, continuità di stato allo staging.

Nessuna soluzione del problema staging-in-tangente-lineare qui (Step 8).

**Avviso:** ciclo più grande del progetto finora (due fasi, prima
modifica di codice testato, primo assemblaggio end-to-end). Parte più
incerta (raggiungibilità fisica) già de-rischiata sopra; possibili round
aggiuntivi se emergono altre sorprese — verranno segnalati esplicitamente.

#### 3bis. Addendum (reviewer) — Fase B da rivedere, 5 problemi concreti
risolti prima di passare al coder

Verdetto reviewer: **da rivedere**. Punti e risoluzioni:

**1. Isp/spinta da vuoto — rischio di "aggiustare finché torna".**
Riconosciuto: la scelta iniziale (vuoto perché il totale cade in range,
SL perché no) era presentata come ovvia senza mostrare perché. Risolto
mantenendo Isp/spinta da vuoto ma con motivazione esplicita e onesta,
non un post-hoc: riportare ESPLICITAMENTE ANCHE il numero con Isp SL
(≈8764 m/s) nel report finale come limite inferiore di sensibilità, e
giustificare la scelta del vuoto con un argomento indipendente dal
risultato (non "perché torna"): lo Stadio 1 con questi dati (spinta
vuoto, kick standard) brucia fino a h≈79-80 km (vedi punto 2 sotto),
quota alla quale la pressione atmosferica e' gia' una frazione minima
di quella al livello del mare per la maggior parte della durata del
bruciamento — le prestazioni vuoto sono percio' una media
pesata-nel-tempo piu' rappresentativa della SL per un modello a Isp
costante. Riportare entrambi i numeri (SL e vuoto) rende la scelta
trasparente invece di nascosta.

**2. Payload — verificare che tutti i numeri vengano dalla STESSA riga
della fonte (variante espendibile, non riutilizzabile).** Task esplicito
per il coder: citare nel codice, per ciascun numero (m_prop, m_strut,
spinta, Isp di entrambi gli stadi, payload), che provengono dalla stessa
tabella "Falcon 9 Block 5, configurazione espendibile" della fonte
Wikipedia — non assumere, verificare e citare esplicitamente.

**3. Identità di chiusura delle perdite — derivata e verificata,
DA DIMOSTRARE anche nel codice.** Il reviewer ha derivato (io l'ho
riverificata algebricamente, stessa identità):

```
d(v)/dt = a(t)*cos(theta(t)-gamma(t)) - g(h)*sin(gamma(t))
```

dove `v=sqrt(vx²+vh²)`, `gamma=atan2(vh,vx)` — i termini centripeti si
CANCELLANO esattamente nel modulo della velocità (non compaiono),
quindi `∫g(h)·sin(gamma)dt` come "perdita gravitazionale" resta valido
anche con la correzione Fase A, per ENTRAMBE le fasi. Ma va dimostrato
nel codice/test (confronto per differenze finite di v(t) simulato
contro l'identità sopra, non solo riusato per analogia dal gravity
turn), come richiesto dal reviewer.

**Punto più serio: il budget di perdite è STRETTO, verificato che lascia
poco margine.** `Δv_ideale (9138.15) - v_orbitale (7788.49) = 1349.66
m/s` di budget totale perdite (gravità+drag+manovra, tutti ≥0, sommano
ESATTAMENTE a questo). Il range di letteratura citato in CLAUDE.md
(1.0-1.5 gravità + 0.1-0.4 drag = 1.1-1.9 totale) ha l'estremo massimo
(1.9) che ECCEDE il budget disponibile (1.35) — impossibile per
costruzione se preso alla lettera come somma di intervalli indipendenti.
**Task obbligatorio per il coder, PRIMA di interpretare i numeri contro
la letteratura:** calcolare i quattro termini (Δv ideale, perdita
gravità, perdita drag, perdita manovra — quest'ultima con formula
diretta `∫a(t)·(1-cos(theta-gamma))dt`, MAI per differenza/residuo) e
verificare che sommino esattamente a `v_finale - v_iniziale` entro la
precisione dell'integratore. Solo dopo, confrontare i singoli termini
con la letteratura — se gravità+drag simulati si avvicinano o superano
1.35 km/s, e' un segnale da investigare esplicitamente (possibile
contributo dell'asimmetria nota del gravity turn, Fase A punto
"Asimmetria dichiarata"), non da far sparire nell'etichetta "perdita di
manovra".

**4. Continuità Stadio1→Stadio2 — task espliciti mancanti, aggiunti
ora:** `vx0_stadio2 = v_fine_stadio1*cos(gamma_fine_stadio1)`,
`vh0_stadio2 = v_fine_stadio1*sin(gamma_fine_stadio1)`; massa di
partenza Stadio 2 = massa EFFETTIVA a fine Stadio 1 (dallo stato di
`solve_ivp`, non nominale) MENO `m_strut_1` — stesso principio già
stabilito per lo staging in fase atmosferica (Step 5, `staging.py`).
Test dedicato che verifica separatamente: (a) la conversione v,gamma→
vx,vh; (b) la sottrazione di m_strut_1; non solo "nessun NaN".

**5. Guess iniziale — l'etichetta "verificato" era INGANNEVOLE, corretto
ora con verifica sulle equazioni COMPLETE.** Il guess A0≈0.22/B0≈-0.0004
era stato verificato PRIMA che il reviewer scoprisse il termine mancante
in dvx/dt (con la correzione "a metà"). **Riverificato da me con le
equazioni complete ora in `guida_esoatmosferica.py`: quel guess NON
converge** (`ier=5`, residuo 33.4). Ho fatto una nuova ricerca numerica
(griglia larga + rifinimento locale): il target esatto (vx=7788.49,
vh=0) **non risulta raggiungibile entro la tolleranza stretta di
`risolvi_coefficienti_tangente_lineare`** con questo `tf` fisso — il
punto più vicino trovato ha uno scarto residuo di **~33-36 m/s
(~0.4-0.5% del target)**, non zero, a `A≈-0.04, B≈0.0009` circa (regione
di best-fit, non un punto esatto).

**Decisione di design (nuova, necessaria):** per Fase B NON si userà
`risolvi_coefficienti_tangente_lineare` con la sua tolleranza stretta
(1e-6, pensata per il caso astratto di Step 6 dove convergenza esatta
era verificata) — si userà invece una ricerca ai minimi quadrati
(`scipy.optimize.minimize`, stesso principio già usato da me per
l'indagine pre-piano) che trova il punto più vicino raggiungibile, e si
**riporta onestamente lo scarto residuo come risultato**, non si forza
una convergenza artificiale. Questo è coerente con lo spirito CLAUDE.md:
un piccolo scarto residuo (~0.4%) con un `tf` fisso e solo 2 gradi di
libertà (A,B) per centrare un bersaglio a 2 coordinate mentre l'altitudine
finale resta libera è un limite onesto e atteso del modello, non un
errore da nascondere — va dichiarato esplicitamente nel report finale
del ciclo, non forzato a zero.

**Prossimo:** coder implementa Fase B secondo questo design consolidato.

#### 4. Esito Fase A (reviewer + coder + critic-ingegnere)

- **reviewer**: verdetto iniziale **da rivedere** — ha derivato le
  equazioni polari complete e dimostrato che la mia proposta iniziale
  (correggere solo dvh/dt) violava conservazione di energia e momento
  angolare (mancava il termine gemello -vh·vx/(R+h) in dvx/dt). Ho
  verificato numericamente (Bash): correzione completa conserva E/L a
  precisione di macchina (~1e-14), quella a metà li viola dell'1-2% in
  200s. Design corretto integrato nell'addendum 1bis prima del coder.
- **coder**: equazioni corrette implementate esattamente. Tutti e 10 i
  test esistenti aggiornati (non solo 4, come inizialmente sottostimato)
  + 2 nuovi test di conservazione. 50/50 verdi. Ha segnalato onestamente
  una conseguenza non prevista dal piano (ValueError da dominio invalido
  durante l'esplorazione di fsolve, dato che g(h) ora è ricalcolata ad
  ogni passo invece di essere congelata) e l'ha risolta con un fallback
  esplicito e documentato (residuo costante finito invece di crash).
  Riverificato che il fenomeno delle radici multiple di Step 6 si
  ripresenta con la nuova fisica (non assunto per analogia).
- **critic-ingegnere**: ricalcolo indipendente dei check di conservazione
  con condizioni iniziali proprie (2 metodi separati, ~1e-11-1e-15 di
  variazione relativa, confermato). Nessun altro modulo toccato
  (verificato via git diff --stat). test_2 confermato non tautologico.
  Ha stress-testato il fallback ValueError→residuo costante con 30 guess
  diversi: 0 convergenze silenziose sbagliate, ma dimostrato che il
  residuo costante NON fornisce gradiente utile a fsolve dentro la zona
  invalida (è un muro, non un pendio — la vera rete di sicurezza è il
  controllo combinato ier+residuo+doppio guess già presente). Nota
  esplicita: il check delta-v (9138.15 m/s calcolato a mano) NON è
  ancora verificabile dal codice reale — da confermare in Fase B, non
  assunto valido solo perché coerente a mano. 50/50 test confermati
  indipendentemente. Nessuna violazione.
- **Correzione docstring (io):** la nota sul fallback ValueError
  suggeriva che "guidasse" l'ottimizzatore — corretto per riflettere la
  scoperta del critic-ingegnere (muro, non pendio). 50/50 confermati
  dopo la modifica (solo docstring, nessuna modifica di logica).

**Fase A: COMPLETATA.**

#### 5. Esito Fase B (coder + critic-ingegnere, pipeline snellita)

- **coder**: implementato `lanciatore/validazione.py` + `tests/test_validazione.py`
  secondo il design consolidato (addendum 3bis). Dati Falcon 9 Block 5
  (espendibile, Wikipedia) citati per ciascun numero. Continuità
  Stadio1→Stadio2 con funzioni isolate e testate separatamente
  (conversione v,γ→vx,vh; sottrazione massa strutturale dalla massa
  EFFETTIVA). Stadio 2 risolto con `scipy.optimize.minimize` (non il
  root-finder a tolleranza stretta di Step 6, motivato esplicitamente:
  il target esatto non è raggiungibile con `tf` fisso). Identità di
  chiusura delle perdite implementata con ri-integrazione densa
  indipendente (scoperta non prevista dal piano: i punti radi
  dell'integrazione a passo adattivo non bastavano per una quadratura
  accurata — risolto senza reimplementare fisica). 65/65 test verdi.
- **Numeri chiave** (tutti calcolati in codice, mai a mano):
  Δv ideale totale (Isp vuoto) = **9138.15 m/s** (dentro 9.1-10.0 km/s);
  Δv ideale sensibilità (Isp1 SL) = 8763.81 m/s (sotto range, riportato
  esplicitamente, non nascosto); v orbitale target (200 km) = 7788.49
  m/s; v finale raggiunta = 7755.09 m/s (scarto 34.4 m/s, 0.44%,
  riportato onestamente, non forzato a zero); perdita gravità = 1271.13
  m/s; perdita drag = 16.73 m/s; perdita manovra = 95.21 m/s; chiusura
  identità: scarto ~1.5e-4 m/s (~2e-8 relativo).
- **critic-ingegnere**: ricalcolo INDIPENDENTE con metodi diversi
  (integratore DOP853 invece di RK45, quadratura Simpson invece di
  trapezoide, script separati che non riusano il codice del progetto) —
  tutti i numeri coincidono esattamente. Identità di chiusura
  soddisfatta a precisione molto più stretta della tolleranza di test.
  Continuità Stadio1→Stadio2 verificata con test non tautologici (valori
  sintetici indipendenti). Nessun altro modulo toccato oltre a
  `guida_esoatmosferica.py` (Fase A). 65/65 confermati indipendentemente.
  **Nessuna violazione**, ma due punti di attenzione segnalati esplicitamente
  (non bloccanti, da tenere presente):
  1. **Margine stretto sul check delta-v**: 9138.15 m/s è dentro il range
     ma solo 38 m/s (0.42%) sopra il limite inferiore (9100 m/s) — il
     risultato dipende in modo sensibile dalla scelta Isp-vuoto vs
     Isp-SL. La validazione passa, ma per un margine stretto, non
     ampio — da tenere presente nei confronti futuri (Step 8/9), non da
     presentare come un fit comodo.
     **Nota dell'utente (2026-08-16, aggiunta anche in CLAUDE.md come
     disclaimer permanente):** il confronto col benchmark 9.1-10.0 km/s
     è un controllo di **plausibilità in un range realistico generico**,
     NON una riproduzione precisa del Falcon 9 reale. Il modello esclude
     per costruzione il bonus di rotazione terrestre (dichiarato fuori
     scope da CLAUDE.md fin dall'inizio) — per un lancio verso est da
     Cape Canaveral (28.5°N) questo bonus vale **≈409 m/s**
     (`465.1 m/s * cos(28.5°)`), calcolato esplicitamente, **maggiore
     del margine di 38 m/s con cui il check è stato superato**. Il check
     resta valido come verifica di ordine di grandezza, non come
     confronto quantitativo diretto con le prestazioni di missione di
     un lancio reale specifico da un sito equatoriale/subtropicale.
  2. **Drag basso (16.7 m/s) plausibile ma al limite del range di
     letteratura (100-400 m/s)**: il critic-ingegnere ha verificato
     indipendentemente il profilo di pressione dinamica dello Stadio 1
     — Max-Q ≈32 kPa a t≈55.7s, h≈10.3km, v≈464 m/s, **sorprendentemente
     vicino ai dati pubblici reali del Falcon 9** (Max-Q reale ≈30-35
     kPa) — buona convalida indiretta del modello drag/atmosfera
     proprio nella banda 0-25km dove il fit di Curtis è più accurato
     (Step 1). Circa il 33% della durata della Fase B (gravity turn)
     avviene sopra i 30 km, dove la densità è già minima. Precedente di
     letteratura citato dal critic-ingegnere: grandi booster ad alto T/W
     (es. Saturn V, drag loss riportata ~40 m/s) hanno perdite da drag
     anch'esse sotto il range 100-400 m/s dichiarato in CLAUDE.md —
     quel range e' probabilmente tarato su profili di ascesa meno
     aggressivi (T/W più basso). **Conclusione: 16.7 m/s è fisicamente
     plausibile per un veicolo ad alto T/W come questo, non un errore**,
     ma resta un caso limite del range dichiarato — investigato,
     spiegato, non nascosto, coerente con la regola CLAUDE.md sui
     risultati inattesi.
  3. Nota di processo dello Step 1 non ancora chiusa: la verifica
     esplicita dell'effetto del range di validità dell'atmosfera
     (0-25km, errore crescente sopra i 30km) sul delta-v da drag era
     richiesta come "da verificare esplicitamente nello Step 6" nella
     docstring originale di `atmosfera.py` (numerazione step precedente
     alla riorganizzazione) — l'analisi Max-Q di cui sopra la soddisfa
     nella sostanza ma andrebbe formalizzata in un punto più visibile
     (es. Step 10, validazione e limiti) invece di restare solo in
     questo log.

**Step 7: COMPLETATO.** Il "check di validazione obbligatorio" di
CLAUDE.md (delta-v vs benchmark 9.1-10.0 km/s) è ora eseguibile dal
codice reale, non solo calcolato a mano — e passa, con margine stretto
ma reale. Questo ciclo ha anche prodotto la prima modifica di codice già
testato nel progetto (termine centripeto, Fase A), scoperta e corretta
con lo stesso rigore (doppia verifica indipendente, derivazione completa,
dimostrazione di conservazione) usato per ogni scoperta precedente.
Prossimo step proposto: Step 8 (estensione — staging durante la fase di
guida a tangente lineare, che riuserà il safeguard a doppio guess già
costruito qui, come annotato in memoria).

---

### 2026-08-16 — Ciclo 8 (piano scritto direttamente dall'orchestratore,
indagine numerica pre-piano)

**Nota di processo:** pipeline snellita confermata: scrivo io il piano,
salto `planner`/`optimizer`/`reporter`, mantengo `reviewer` e
`critic-ingegnere`.

**Step:** 8 — Estensione: staging durante la fase di guida a tangente
lineare.

**Promemoria onorato (nota utente + memoria persistente, Ciclo 6):** il
problema inverso A/B ha radici multiple/spurie per un singolo segmento;
atteso che il fenomeno si ripresenti con più segmenti concatenati. Il
safeguard a doppio guess già in `risolvi_coefficienti_tangente_lineare`
va riusato, non reinventato.

#### 1. Design (verificato numericamente prima di scrivere il piano)

**Scenario dimostrativo:** un segmento di guida a tangente lineare
pianifica A1,B1 per un target finale assumendo (erroneamente) di
bruciare per l'intera durata nominale — ma lo staging avviene PRIMA del
previsto, invalidando quel piano. Il nuovo segmento deve ricalcolare
A2,B2 dallo stato EFFETTIVO di handoff per raggiungere comunque lo
stesso target. Scelto invece di un target intermedio arbitrario perché
DIMOSTRA (non solo asserisce) il motivo per cui questo step esiste.

**Parametri verificati (Bash, convergenza robusta):**
- Segmento 1: m0=1000 kg, spinta=40000 N, isp=300 s, h0=5000 m (margine
  quota positivo, stesso principio di Fase A/Ciclo 7). Target finale:
  vx=1800 m/s, vh=0 (entro capacità ideale Tsiolkovsky ≈2309 m/s del
  solo segmento 1). tf1 nominale=40s. `A1=2.19215, B1=-0.081227`
  (convergenza identica da 5/6 guess diversi).
- **Staging anticipato a t=15s**: stato allo staging `x=2412.78,
  h=8077.61, vx=369.26, vh=408.74, m=796.057`. `m_strut1=100 kg`
  (provvisorio) → `m_ignizione2=696.057 kg`.
- **Prova che A1,B1 non sono più validi:** integrando il segmento 2 con
  A1,B1 invariati sotto i nuovi parametri, si ottiene `vx=771.37,
  vh=565.59` — lontano dal target (1800, 0).
- Segmento 2: spinta=40000 N, isp=320 s, tf2=25s. Ri-risolto dallo
  stato EFFETTIVO di handoff: `A2=1.859568, B2=-0.145024`, target
  raggiunto ESATTAMENTE (vx=1800.000, vh=0.00000).
- **Fenomeno radici multiple RICONFERMATO** (non solo atteso per
  analogia): guess `(A0=0.5, B0=-0.01)` converge a una radice DIVERSA
  ma ugualmente valida (`A2=-2.29081, B2=0.149645`, stesso target
  raggiunto esattamente).

#### 2. Nuovo modulo

`lanciatore/staging_esoatmosferico.py`: (1) risolve A1,B1 nominali; (2)
integra segmento 1 solo fino al tempo di staging; (3) applica
discontinuità di massa (stato effettivo meno massa strutturale, stesso
principio di `staging.py`); (4) ri-risolve A2,B2 dallo stato di
handoff; (5) ritorna entrambi i segmenti + i coefficienti di entrambi i
piani (originale e ricalcolato) per il confronto esplicito. Riuso
totale di `guida_esoatmosferica.py`, nessuna modifica a moduli
esistenti.

#### 3. Criteri di verifica

1. Continuità di stato allo staging (x,h,vx,vh esatti; massa che
   differisce esattamente di `m_strut1`).
2. **Dimostrazione esplicita** che A1,B1 non validi post-staging
   (errore >10% integrando con i vecchi coefficienti) PRIMA di
   verificare che il ricalcolo funzioni.
3. A2,B2 ricalcolati raggiungono il target entro la tolleranza del
   risolutore (riusata, non allargata).
4. Documentazione del fenomeno radici multiple in questo contesto
   multi-segmento.
5. Nessun NaN/inf sulla traiettoria concatenata.

Nessun dato reale in questo step (scenario di test, come Step 5).
Nessun nuovo check delta-v (resta quello di Step 7).

**Prossimo:** passare al sub-agente reviewer per il controllo critico
di questo piano, poi al coder.

#### 4. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **da rivedere**, tre problemi concreti, tutti
risolti e verificati numericamente prima del coder:

**1. Lo staging deve essere un vero evento ODE terminale, non un
cutoff scriptato.** Riformulato lo scenario: il Segmento 1 non ha
propellente per l'intera durata nominale ipotizzata in fase di
pianificazione (40s) — porta realisticamente solo
`m_prop1_reale = mdot1*15s ≈ 203.94 kg` (mdot1≈13.596 kg/s), quindi
esaurisce il propellente a t≈15s per un vincolo fisico reale del
veicolo, non per un taglio arbitrario. Soglia di massa
`m_soglia1 = m0 - m_prop1_reale ≈ 796.06 kg`. Implementato come evento
terminale genuino (`y[4]-m_soglia1`, `terminal=True`, `direction=-1`),
stesso pattern di `evento_fine_propellente_2d` in `staging.py`, con
`solve_ivp` chiamato DIRETTAMENTE (bypassando `integra_tangente_lineare`,
che non espone `events`) — stesso approccio già usato in `staging.py`
per gli stadi 2..N.

**2. Regola di seeding deterministica per il guess del Segmento 2:**
usare direttamente `A0=A1, B0=B1` (i coefficienti convergenti del
segmento precedente) come guess per ri-risolvere il segmento
successivo — verificato numericamente che converge in modo affidabile
(3/3 prove ripetute) alla radice fisicamente continua
(`A2=1.85957, B2=-0.145024`), non a quella spuria. Coerente col
principio già scritto nella docstring di
`risolvi_coefficienti_tangente_lineare` ("PEG non riparte mai da un
default generico, riusa la soluzione del ciclo precedente"). Nuovo
criterio di verifica esplicito: la funzione di staging usa SEMPRE
questa regola (non un guess arbitrario), e un test verifica la stabilità
ripetendo la risoluzione 3 volte confermando convergenza identica.

**3. Metrica di errore e continuità temporale precisate per la "prova
di invalidità" (criterio 2).** Continuare il piano ORIGINALE senza
ricalcolo significa continuare il MEDESIMO profilo `theta(t)` con
tempo assoluto continuo (non un orologio che riparte da zero al
segmento 2) — equivalente a integrare il segmento 2 con
`A0_continuato = A1 + B1*t_stage`, `B0_continuato = B1` (stesso B,
A traslato dell'offset temporale). **Ricalcolato con questa convenzione
corretta** (la prima stima nel piano, `vx=771/vh=566`, era gonfiata da
un artefatto di reset del tempo, non solo dalla discontinuità fisica):
risultato onesto **vx=2045.76, vh=-41.90** contro target (1800, 0) —
errore euclideo 249.3 m/s, **errore relativo 13.85%** (sopra la soglia
10%, criterio ancora soddisfatto ma con un numero corretto, non
gonfiato). Metrica precisata: errore euclideo in (vx,vh) normalizzato
sulla norma del target (non errore per-componente, indefinito per
`vh_target=0`).

**Confermato senza modifiche (verificato positivamente dal reviewer):**
continuità Stadio1→Stadio2 è un riporto diretto (x,h,vx,vh, nessuna
conversione trigonometrica necessaria, a differenza dello Step 7 dove
si passava da stato polare v,γ a cartesiano vx,vh) — caso più semplice
di `staging.py`, da non confondere aggiungendo l'helper di conversione
dello Step 7 dove non serve.

**Prossimo:** coder implementa secondo questo design consolidato.

#### 5. Esito ciclo (coder + critic-ingegnere, pipeline snellita)

- **coder**: implementato `lanciatore/staging_esoatmosferico.py` +
  `tests/test_staging_esoatmosferico.py` esattamente secondo il design
  consolidato. Nessuna ambiguità residua (piano+addendum sufficientemente
  precisi). 73/73 test verdi (65 pregressi + 8 nuovi). Tutti i numeri
  coincidono con la verifica pre-piano dell'orchestratore: A1=2.1921496,
  B1=-0.0812265; evento genuino a t_stage=15.0s esatto (status=1, non
  t_span esaurito); A2=1.8595704, B2=-0.1450243 (target raggiunto,
  vx=1799.9999996, vh≈-1.8e-8); errore "senza ricalcolo" 13.85%; radice
  spuria confermata (A=-2.290804, B=0.149645, target comunque raggiunto).
- **critic-ingegnere**: ricalcolo INDIPENDENTE di ogni numero chiave
  (script separato) — coincidenza esatta. Vincolo "Multistadio: eventi
  di staging come eventi terminali ODE" verificato esplicitamente
  (status==1, t_events=[15.0], non un cutoff scriptato). Regola di
  seeding A0=A1/B0=B1 verificata nel codice e ripetuta 5 volte
  (convergenza identica bit-per-bit). Fenomeno radici multiple
  riconfermato con un TERZO guess indipendente non presente nel piano
  (-1.0, 0.05), anch'esso convergente alla stessa radice spuria — non
  un artefatto di un guess specifico. Nessun modulo esistente
  modificato (git diff vuoto). 73/73 confermati indipendentemente.
  Nessuna violazione. Ha ri-segnalato (non nuovo, pre-esistente e già
  documentato) il margine stretto del delta-v di missione dello Step 7.
- **Pulizia (io, non optimizer):** test già ben organizzati, nessuna
  correzione necessaria.

**Step 8: COMPLETATO.** Il fenomeno delle radici multiple scoperto allo
Step 6 si è confermato riemergere in un contesto multi-segmento, esattamente
come anticipato dalla nota utente/memoria persistente — e il safeguard
già costruito (regola di seeding deterministica, non un default
arbitrario) si è dimostrato sufficiente senza dover reinventare nulla.
Anche il reviewer, in questo ciclo, ha trovato un problema reale prima
dell'implementazione: la prima versione del "test di invalidità"
avrebbe gonfiato l'errore per un artefatto di reset del tempo (report
di un ~65% invece del 13.85% reale) — corretto prima di scrivere codice.
Prossimo step proposto: Step 9 (visualizzazione, traiettoria numerica +
animazione).

---

### 2026-08-16 — Ciclo 9 (piano scritto direttamente dall'orchestratore)

**Nota di processo:** pipeline snellita confermata: scrivo io il piano,
salto `planner`/`optimizer`/`reporter`, mantengo `reviewer` e
`critic-ingegnere`. Natura diversa dai cicli precedenti: nessuna nuova
fisica, solo visualizzazione di dati già calcolati e verificati.

**Step:** 9 — Visualizzazione (traiettoria numerica + animazione).

**Dataset scelto:** la traiettoria completa dello Step 7
(`lanciatore.validazione.esegui_validazione()`/`assembla_traiettoria()`)
— il risultato più completo del progetto (Falcon 9, gravity turn +
tangente lineare fino a orbita). Riuso totale, nessuna nuova
integrazione.

#### 1. Design

Nuovo modulo `lanciatore/visualizzazione.py` (nessuna modifica a moduli
esistenti):
1. `estrai_serie_temporali(...)`: converte i risultati grezzi di
   `assembla_traiettoria()` in serie temporali concatenate a TEMPO
   ASSOLUTO CONTINUO attraverso le 3 fasi (Fase A verticale, Fase B
   gravity turn, tangente lineare — sommando gli offset, non
   resettando l'orologio, stesso principio corretto dal reviewer allo
   Step 8), con marcatori espliciti delle transizioni di fase.
2. `grafico_quota_velocita_tempo(...)`: h(t), v(t) modulo, m(t) —
   sottografici con linee verticali alle transizioni, salvato PNG.
3. `grafico_traiettoria(...)`: piano (x,h), salvato PNG.
4. `anima_traiettoria(...)`: `FuncAnimation` su (x,h), salvata GIF via
   `PillowWriter` (Pillow già installato, nessuna nuova dipendenza).
5. Riepilogo numerico testuale (delta-v ideale/raggiunto, scomposizione
   perdite — riusa `calcola_delta_v_ideale`/`scompone_perdite`, GIA'
   calcolati e verificati allo Step 7, non ricalcolati qui).

Output in `output/` (nuova cartella, aggiunta a `.gitignore` — artefatti
rigenerabili, stesso trattamento di `venv/`).

#### 2. Criteri di verifica

1. Serie temporali estratte: stessa lunghezza dei dati sorgente, tempo
   monotono crescente attraverso le transizioni (continuità, non reset).
2. Marcatori di transizione coerenti con gli eventi reali (es. fine
   Fase A = evento di kick).
3. Funzioni di grafico/animazione girano senza eccezioni, file salvati
   non vuoti.
4. **Ispezione visiva effettiva** (io, non solo test automatici): genero
   i grafici, li guardo (Read su PNG) prima di dichiarare concluso,
   li invio all'utente (SendUserFile).

**Prossimo:** passare al sub-agente reviewer per il controllo critico
di questo piano, poi al coder.

#### 3. Addendum (reviewer) — recepito prima di passare al coder

Verdetto reviewer: **da rivedere** (non problemi di fisica, ma
ambiguità concrete di estrazione dati che rischiano un bug silenzioso
in un grafico "plausibile a occhio ma quantitativamente sbagliato").
Verificato nel codice riga per riga. Risolto:

**1. Formula esplicita per l'offset di tempo** (i tre `OdeResult` sono
indipendenti, ciascuno riparte da `t=0`):
```
offset_B = fase_verticale.t[-1]
offset_tangente = offset_B + gravity_turn.t[-1]
```

**2. Tabella esplicita di mapping stato→grandezza, DIVERSA per ogni
fase (non un indice uniforme):**

| Fase | Stato | h | v (modulo) | m | x |
|---|---|---|---|---|---|
| Fase A (verticale) | `[h,v,m]` | idx0 | idx1 (diretto) | idx2 | **non esiste, sintetizzare array di zeri** |
| Fase B (gravity turn) | `[x,h,v,gamma,m]` | idx1 | idx2 (diretto) | idx4 | idx0 |
| Tangente lineare | `[x,h,vx,vh,m]` | idx1 | **`sqrt(vx**2+vh**2)`, idx2/idx3** | idx4 | idx0 |

**3. `x` per la Fase A non esiste nei dati grezzi:** sintetizzare
esplicitamente `np.zeros_like(fase_verticale.t)` (ascesa verticale pura,
x=0 per costruzione) — NON un `IndexError`/riuso implicito di un'altra
colonna.

**4. Nuovo criterio di verifica quantitativo (chiusura con i numeri già
verificati allo Step 7, non solo ispezione visiva):** il valore finale
della serie concatenata di velocità-modulo, quota, e massa deve
coincidere ENTRO TOLLERANZA NUMERICA con
`assemblaggio['v_finale']`/`h` finale di `risultato_2`/`m` finale di
`risultato_2` (già calcolati e verificati indipendentemente dal
critic-ingegnere allo Step 7) — e la massa iniziale di ciascun segmento
deve coincidere con `m1_effettiva`/`m_ignizione_2_effettiva`. Non un
nuovo delta-v da validare, ma una chiusura numerica tra dati grezzi già
validati e serie derivata per il plotting.

**5. Discontinuità di massa allo staging (Stadio1→Stadio2): comportamento
ATTESO, non un bug da "correggere".** Concatenando gli array grezzi così
come sono, la serie `m(t)` mostra correttamente un salto verso il basso
nel punto di transizione (nessuna interpolazione). Dichiarato
esplicitamente per il coder: il salto va PRESERVATO nel grafico, non
smussato.

**Prossimo:** coder implementa secondo questo design consolidato.

#### 4. Esito ciclo (coder + critic-ingegnere, pipeline snellita)

- **coder**: implementato `lanciatore/visualizzazione.py` +
  `tests/test_visualizzazione.py` secondo il design consolidato.
  `output/` aggiunto a `.gitignore`. 86/86 test verdi (73 pregressi + 13
  nuovi). 3 file generati: `output/quota_velocita_tempo.png`,
  `output/traiettoria.png`, `output/traiettoria.gif`. Segnalata
  un'imprecisione di formulazione nell'addendum (punto 4, "massa
  iniziale di ciascun segmento" — `m1_effettiva` è in realtà la massa
  FINALE del gravity turn) risolta con l'unica lettura sensata
  (verifica dei due lati della discontinuità), non ignorata.
- **Ispezione visiva (io):** grafici corretti — quota/velocità/massa
  monotone come atteso nella Fase A/B, salto di massa allo staging
  visibile come gradino netto, traiettoria (x,h) con curva del gravity
  turn che si raddrizza dopo lo staging. **Notato un dettaglio nel
  grafico h(t) non ovvio dai soli numeri aggregati dello Step 7: la
  quota sale a ~170km poi ridiscende a ~148km entro fine bruciamento
  del Segmento 2** — verificato non essere un bug di visualizzazione
  (dati confermati identici a `validazione.py`, non toccato).
- **critic-ingegnere**: tutti i criteri dell'addendum verificati con
  ricalcolo indipendente (mapping per fase, offset tempo, salto di
  massa esatto -25600kg, chiusura quantitativa esatta). 86/86
  confermati. **Ha approfondito la scoperta della quota discendente
  fino alle sue conseguenze fisiche**: propagando lo stato finale dello
  Step 7 con la meccanica orbitale standard (nessuna nuova fisica),
  l'orbita osculatrice ha perigeo calcolato a **-62.5 km** (sotto la
  superficie terrestre) — **NON un'orbita LEO stabile**, una traiettoria
  che rientrerebbe prima di completare un giro. **Confermato
  indipendentemente anche da me** (Bash, leggi orbitali standard:
  E=-31 073 249 J/kg, semiasse maggiore=6413.9km, eccentricità=0.01643,
  perigeo=-62.50km, apogeo=148.27km — coincide esattamente col calcolo
  del critic-ingegnere). Causa: il risolutore dello Step 7
  (`risolvi_stadio_2_minimi_quadrati`) vincola solo velocità finale, non
  quota — scelta di scope esplicita di quello step ("l'altitudine
  finale resta libera"), che però ha una conseguenza più seria di quanto
  notato allora: il check delta-v (scalare aggregato) non intercetta
  un'orbita non valida.
- **Pulizia (io, non optimizer):** codice già ben strutturato, nessuna
  correzione necessaria.

**Step 9: COMPLETATO nello scope dichiarato** (grafici e animazione
fedeli ai dati, 86/86 test, output generato e ispezionato). **Ha però
prodotto la scoperta più importante del progetto sul piano dei
risultati** (non solo del processo): la validazione Falcon 9 dello Step
7, per quanto passi il check delta-v obbligatorio, non produce
un'orbita LEO effettivamente stabile. Decisione su come/quando
affrontarlo lasciata all'utente (vedi nota nella roadmap sopra) —
possibilità: (a) risolvere ora, tornando allo Step 7 per vincolare
anche la quota nel solutore (probabilmente serve un terzo grado di
libertà, es. `tf` libero); (b) documentarlo come limite noto nello
Step 10 (già dedicato a "validazione e documentazione dei limiti") e
procedere; (c) altro. Prossimo step proposto: da confermare con
l'utente alla luce di questa scoperta.

---

### 2026-08-17 — Ciclo 10 (correzione mirata a Step 7, scritta
direttamente dall'orchestratore, indagine numerica pre-piano)

**Decisione dell'utente:** tornare allo Step 7 e correggere il
risolutore dello Stadio 2 perché vincoli **quota finale = target** e
**componente verticale della velocità finale = 0**, non più magnitudine/
componente orizzontale della velocità. Verificare poi che lo stato
finale dia un'orbita reale (perigeo sopra la superficie) con lo stesso
calcolo di meccanica orbitale già usato per la diagnosi. Riverificare
Step 8 e 9 dopo la correzione.

#### 1. Indagine numerica pre-piano (fatta prima di scrivere il piano,
stesso principio già seguito per la scoperta del termine centripeto)

**Step 8 confermato indipendente da `validazione.py`:** verificato via
grep degli import di `staging_esoatmosferico.py`/`test_staging_esoatmosferico.py`
— nessun riferimento a `validazione`. Scenario autocontenuto (m0=1000,
parametri di test propri). Non necessita modifiche, solo riconferma
finale della suite.

**Nuovo obiettivo testato numericamente:** `risolvi_stadio_2_minimi_quadrati`
deve minimizzare `(h(tf)-h_target)² + vh(tf)²` invece di
`(vx(tf)-vx_target)² + (vh(tf)-vh_target)²`. Testato con lo stato di
handoff reale (identico a prima, Stadio 1/gravity turn non tocco):
`x1=113793.17, h1=78943.84, vx1=2726.90, vh1=1169.41, m_ignizione2=119370.0`.

**Risultato (guess `A0=0.3, B0=0.0`, converge in modo pulito — altri
guess provati falliscono per `ValueError` di dominio durante
l'esplorazione, comportamento già noto e gestito dal fallback
residuo-costante di Ciclo 7):** `A=0.24410, B=-0.000522`,
`h_finale=200.000 km` (esatto), `vh_finale≈0`, **`vx_finale=7721.04 m/s`**
(contro `vx_orbitale=7788.49 m/s` — 67.4 m/s in meno, il vincolo di
quota "consuma" margine di velocità che prima andava tutto su vx).

**Verifica orbitale (stessa meccanica standard già usata per la
diagnosi, nessuna nuova fisica):** `a=6459.61 km`, `ecc=0.017244`,
**perigeo=-22.78 km, apogeo=200.00 km**. **Migliore della versione
precedente (-62.5 km→-22.8 km, -63% di scarto) ma ANCORA NEGATIVO —
non ancora un'orbita stabile.**

**Esplorazione aggiuntiva (per informare l'utente, non per scegliere
unilateralmente un target diverso):** ripetuto lo stesso calcolo per
`h_target` da 150 a 200 km — il perigeo resta negativo su tutto
l'intervallo (minimo scarto ≈-15 km attorno a 180-190 km, non
sostanzialmente migliore di -22.8 km a 200 km). **Non è un problema di
quale quota target scegliere: è un vero scarto di budget Δv/energia
dello Stadio 2 con questo `tf` fisso e questo stato di handoff** — la
stessa cosa già segnalata come "margine stretto" nel disclaimer
CLAUDE.md del check delta-v (38 m/s sopra il minimo), ora quantificata
in modo più diretto e fisico (mancano di ordine 20-25 km di perigeo,
non solo m/s aggregati).

#### 2. Design della correzione

- `lanciatore/validazione.py`: nuova funzione obiettivo
  `_obiettivo_residuo_quota_verticale` (minimizza `(h(tf)-h_target)² +
  vh(tf)²`), affiancata a quella esistente (NON sostituita — la vecchia
  resta disponibile per confronto/sensibilità, marcata esplicitamente
  come "target precedente, lasciava la quota libera"). `risolvi_stadio_2_minimi_quadrati`
  aggiornata per accettare il nuovo target (h_target invece di
  vx_target) con lo stesso pattern di doppio guess/diagnostica di
  robustezza già presente.
- Nuova funzione `verifica_orbita_stabile(h_finale, vx_finale, vh_finale)`:
  calcola semiasse maggiore, eccentricità, perigeo, apogeo con la
  meccanica orbitale standard (energia specifica + momento angolare,
  MU_TERRA/R_TERRA già citati) — riusabile anche in step futuri (es.
  Step 10). Ritorna esplicitamente se il perigeo è sopra la superficie
  (`perigeo_valido: bool`), senza nascondere il caso negativo.
- `assembla_traiettoria`: aggiornata per usare il nuovo target di
  default, riporta ESPLICITAMENTE nel risultato sia il vecchio residuo
  (vx) sia il nuovo (h), e il risultato di `verifica_orbita_stabile`.
- **Nessuna modifica a `lanciatore/guida_esoatmosferica.py`,
  `lanciatore/guida.py`, `lanciatore/staging_esoatmosferico.py`,
  `lanciatore/staging.py`** — solo `validazione.py`.

#### 3. Criteri di verifica

1. Nuovo target (h,vh) raggiunto entro tolleranza stretta esplicita.
2. `verifica_orbita_stabile` applicata allo stato finale — risultato
   riportato ESPLICITAMENTE (perigeo negativo o positivo, non nascosto).
3. Δv ideale totale (Step 7, invariato — non dipende dal target dello
   Stadio 2) resta in 9.1-10.0 km/s, riconfermato non ricalcolato.
4. **Step 8**: nessuna modifica necessaria (confermato indipendente),
   solo riconferma che la sua suite resta verde.
5. **Step 9**: `visualizzazione.py` consuma `assembla_traiettoria()`
   dinamicamente (nessun valore hardcoded verificato al Ciclo 9) — atteso
   che i test di chiusura quantitativa restino verdi automaticamente
   con i nuovi numeri, MA i 3 file di output (PNG×2 + GIF) vanno
   RIGENERATI (i vecchi mostrano la traiettoria col vecchio target) e
   re-ispezionati.
6. Suite completa verde.

**Onestà del risultato atteso:** il perigeo resta negativo
(-22.8 km) anche con la correzione. Questo va riportato esplicitamente
all'utente come esito, non nascosto né presentato come "risolto" — la
correzione è comunque quella concettualmente giusta (vincolare quota e
verticalità invece di lasciare la quota libera), e riduce lo scarto in
modo sostanziale, ma non lo azzera con questa configurazione di
veicolo/target/tf.

**Prossimo:** reviewer sul design, poi coder.

#### 4. Addendum (reviewer + tre round di indagine numerica con l'utente)
— design finale, recepito prima del coder

**Verdetto reviewer: da rivedere.** Trovato un punto concettuale
importante: con `h_target` fissata e `vh=0`, quel punto è
l'**APOGEO** dell'orbita (perché il semiasse maggiore risultante `a`
è minore del raggio target), non il perigeo — il vincolo scelto non
controlla affatto la grandezza che si vuole validare. Questo spiega
perché lo sweep di `h_target` (150-200km) non cambiava molto il
perigeo. Altri problemi del reviewer: formule di `verifica_orbita_stabile`
da esplicitare con fonte, test di regressione contro i numeri già
noti, evitare dead code sulla vecchia funzione obiettivo, normalizzare
i termini dell'obiettivo (metri² vs (m/s)² altrimenti mal
condizionato), rivedere la tolleranza di chiusura di `scompone_perdite`
con la nuova forma di traiettoria.

**Tre round di indagine numerica aggiuntiva con l'utente (tutti fatti
PRIMA di scrivere il design finale):**

1. **Obiettivo auto-consistente** (velocità circolare sulla quota
   EFFETTIVAMENTE raggiunta, non fissata a priori — l'alternativa
   proposta dal reviewer): testato, converge a perigeo=-15.86km — il
   migliore trovato con target fisso, ma ancora negativo.
2. **`tf` reso libero** (terza incognita, sistema A/B/tf a 3 equazioni
   invece di 2 sovradeterminate) su richiesta esplicita dell'utente:
   testato con `fsolve` e `minimize`, **converge sempre a
   `tf≈tf2_max`** (bruciamento completo, 322.38s) — l'ottimizzatore
   vuole già tutto il propellente disponibile. **`tf` libero NON aiuta:
   conferma che è un vero deficit di propellente/energia, non un
   problema di formulazione del target o di tempistica.**
3. **Sweep del payload** (per completezza, non per scegliere
   unilateralmente): a payload=22800kg (max dichiarato) perigeo=-21km;
   a payload=20000kg (-12%) perigeo=+199km (quasi circolare); resta
   valido fino a ~7500kg, degrada di nuovo sotto 5000kg.
4. **Bonus di rotazione terrestre** (409 m/s, già quantificato nel
   disclaimer CLAUDE.md dal Ciclo 7) su richiesta esplicita
   dell'utente: aggiunto post-hoc a `vx_finale` (stesso principio con
   cui si applica realmente — offset costante di velocità inerziale
   dovuto al sito di lancio in rotazione, non richiede rieseguire la
   traiettoria), **il perigeo passa da -21km a +199.9km** — il deficit
   trovato (67.4 m/s in vx) è **interamente coperto e superato** dal
   bonus (408.74 m/s, con margine di 341.3 m/s, motivo per cui
   l'orbita risultante diventa ellittica, apogeo ~1495km, se non si
   ribilancia anche il target).

**Decisione finale dell'utente:** implementare il fix richiesto
originariamente (target quota+vh=0, NON l'obiettivo auto-consistente,
NON `tf` libero — nessuno dei due migliora sostanzialmente il
risultato), **payload al valore massimo dichiarato (22800kg, nessun
aggiustamento)**, e **documentare esplicitamente che il perigeo
residuo (~-21/-23km) è interamente spiegato dal bonus di rotazione
terrestre già escluso per scope fin dall'inizio del progetto** — non
un limite aperto, una conferma quantitativa di coerenza del modello.

**Design finale consolidato per `lanciatore/validazione.py`:**

1. **Nuova coppia di funzioni, la vecchia resta invariata e testata
   separatamente (nessun dead code, reviewer punto 4):**
   `_obiettivo_residuo_quota_verticale(coefficienti, m0, spinta, mdot,
   tf, x0, h0, vx0, vh0, h_target, vh_target)` — obiettivo
   NORMALIZZATO (reviewer punto 6):
   `((h(tf)-h_target)/h_target)**2 + (vh(tf)/v_target_circolare)**2`
   dove `v_target_circolare = sqrt(MU_TERRA/(R_TERRA+h_target))` (per
   dare scala fisicamente sensata al termine di velocità). Fallback
   ValueError→1e12 stesso principio già stabilito.
   `risolvi_stadio_2_target_quota(...)` — stesso pattern di doppio
   guess/diagnostica di `risolvi_stadio_2_minimi_quadrati` (che RESTA
   nel modulo, non toccata, usata da un proprio test standalone per
   restare "raggiungibile").
2. **`verifica_orbita_stabile(h, vx, vh, bonus_rotazione=0.0)`** — nuova
   funzione, riusabile. Formule esplicite (fonte: Curtis, *Orbital
   Mechanics for Engineering Students*, già citato dal progetto dallo
   Step 1 per MU_TERRA/R_TERRA):
   ```
   r = R_TERRA + h
   vx_eff = vx + bonus_rotazione   # bonus applicato SOLO alla componente orizzontale
   L = r * vx_eff                  # momento angolare specifico
   E = (vx_eff**2 + vh**2)/2 - MU_TERRA/r   # energia specifica
   a = -MU_TERRA/(2*E)             # semiasse maggiore
   e = sqrt(1 + 2*E*L**2/MU_TERRA**2)   # eccentricità
   r_perigeo = a*(1-e); r_apogeo = a*(1+e)
   ```
   Ritorna dict con `semiasse_maggiore`, `eccentricita`, `perigeo`
   (quota, m), `apogeo` (quota, m), `perigeo_valido` (bool,
   `perigeo >= 0`). Parametro opzionale `bonus_rotazione` (default 0,
   il caso rigoroso "come dichiarato in CLAUDE.md, niente rotazione
   terrestre"; passare `408.74` per la sensibilità con rotazione).
3. **Test di regressione (reviewer punto 3, oracolo dai numeri già
   due volte validati indipendentemente al Ciclo 9):**
   `verifica_orbita_stabile(148050.2, 7755.084, -8.241)` deve dare
   perigeo≈-62.50km, apogeo≈148.27km (stato vecchio target);
   `verifica_orbita_stabile(200000.0, 7721.045, 0.0)` deve dare
   perigeo≈-22.78km, apogeo≈200.00km (stato nuovo target, valori
   dell'indagine pre-piano dell'orchestratore, guess pulito
   `A0=0.3,B0=0.0`, convergenza a precisione di macchina `fun=2.4e-14`).
4. **`assembla_traiettoria`**: passa a usare
   `risolvi_stadio_2_target_quota` di default. Riporta ESPLICITAMENTE
   nel risultato sia il residuo vecchio (vx) sia il nuovo (h), più
   `verifica_orbita_stabile(...)` SENZA bonus (caso rigoroso, ci si
   aspetta `perigeo_valido=False`) E con bonus (`bonus_rotazione=408.74`,
   ci si aspetta `perigeo_valido=True`) — entrambi nel dict di ritorno,
   nessuno nascosto.
5. **Test esplicito che documenta il limite come "spiegato, non
   aperto"** (reviewer punto 8, stesso principio già usato con
   successo per lo scarto dell'11km in `test_atmosfera.py`, citato da
   CLAUDE.md come esempio corretto): asserire `perigeo_valido == False`
   SENZA bonus sullo scenario Falcon 9 reale, poi `perigeo_valido ==
   True` CON bonus — commento esplicito che collega i due, non solo
   un numero isolato.
6. **`scompone_perdite`**: riverificare esplicitamente (non assumere)
   che la tolleranza di chiusura resti stretta con la nuova traiettoria
   (A,B diversi). Se non regge con la griglia attuale, ritarare i
   punti di valutazione, documentando il nuovo valore.

**Impatto su Step 8:** nessuna modifica di codice (confermato
indipendente da `validazione.py`), solo riconferma della suite.

**Impatto su Step 9:** `visualizzazione.py` consuma `assembla_traiettoria()`
dinamicamente — atteso che i test di chiusura quantitativa restino
verdi automaticamente. I 3 file di output vanno RIGENERATI. **Aspettativa
qualitativa esplicita per la re-ispezione (reviewer):** il nuovo
grafico h(t) NON dovrebbe più mostrare il pattern "sale a 170km poi
scende a 148km" del Ciclo 9 — ora `vh(tf)→0`, quindi la salita
dovrebbe appiattirsi verso il target senza ridiscesa marcata.

**Aggiornamento disclaimer CLAUDE.md:** rafforzare la nota già presente
(Step 7) aggiungendo la chiusura quantitativa: lo scarto di perigeo
(-21/-23km) trovato al Ciclo 10 è interamente spiegato e superato dal
bonus di rotazione terrestre già lì quantificato (409 m/s > 67 m/s di
deficit in velocità).

**Prossimo:** coder implementa questo design consolidato.

#### 5. Esito ciclo (coder + critic-ingegnere, pipeline snellita)

- **coder**: implementato tutto il design consolidato in
  `lanciatore/validazione.py` (nuove funzioni aggiunte, `risolvi_stadio_2_minimi_quadrati`
  esistente NON toccata, byte-identica — verificato dal critic-ingegnere).
  91/91 test verdi (85 pregressi + 6 nuovi). Nuovo target raggiunto a
  precisione di macchina (h_finale=200.000km esatto, vh_finale≈0).
  Perigeo SENZA bonus = -22.78km (`perigeo_valido=False`), CON bonus
  rotazione = +200.00km esatto (`perigeo_valido=True`) — entrambi
  riportati esplicitamente nel risultato, nessuno nascosto. Identità di
  chiusura di `scompone_perdite` riverificata con la nuova traiettoria:
  regge senza bisogno di ritarare la griglia (residuo ~1.5e-4 m/s,
  ~1.9e-8 relativo). Rigenerati i 3 file di output dello Step 9 — il
  pattern "sale poi scende" del Ciclo 9 è sparito (resta un lieve
  overshoot di ~2km, tre ordini di grandezza più piccolo, atteso: il
  target vincola solo lo stato a `t=tf`, non la monotonia nel mezzo).
  Segnalata onestamente un'ambiguità sulla copertura di test della
  vecchia funzione dopo il cambio di default, risolta con un nuovo test
  standalone dedicato invece di lasciarla diventare dead code.
- **critic-ingegnere**: ricalcolo INDIPENDENTE (script separato) di
  entrambi i casi oracolo (perigeo -62.50/-22.77km, coincidenza a
  centesimi di km) e del caso con bonus (perigeo +200.00km esatto,
  atteso per costruzione dato che vh=0 implica un apside a quella
  quota). Confermato che il nuovo target è raggiunto a precisione di
  macchina nell'assemblaggio reale. Confermato `risolvi_stadio_2_minimi_quadrati`
  byte-identica e genuinamente testata (non dead code). Identità di
  chiusura riverificata indipendentemente. Step 8/9 confermati non
  toccati e verdi. Nessun altro modulo toccato. 91/91 confermati.
  Δv ideale totale invariato (9138.15 m/s, dentro range, margine
  invariato). **Un solo elemento del piano non ancora eseguito
  segnalato: l'aggiornamento del disclaimer CLAUDE.md con la chiusura
  quantitativa** — completato subito dopo dall'orchestratore (vedi
  sotto). Segnalata anche una fragilità minore di leggibilità (i guess
  di default di `esegui_validazione` puntavano ai vecchi costanti,
  l'ottimizzatore convergeva comunque correttamente ma il nome non
  corrispondeva all'uso) — corretta dall'orchestratore.
- **Chiusura (io):** aggiornato il disclaimer in CLAUDE.md con la
  chiusura quantitativa del Ciclo 10 (nuova sezione "Chiusura
  quantitativa", 2026-08-18). Corretto il default di
  `esegui_validazione` per usare i guess `*_QUOTA` corretti (nessuna
  modifica di logica, solo chiarezza). 91/91 riconfermati dopo la
  modifica.

**Ciclo 10: COMPLETATO.** La scoperta più significativa del progetto
(Step 9: la validazione Falcon 9 non produceva un'orbita stabile) è ora
chiusa con una spiegazione quantitativa verificata due volte in modo
indipendente: il deficit di perigeo non è un difetto del codice o della
fisica implementata, è la conseguenza diretta e ora misurata di una
scelta di scope dichiarata fin dall'inizio del progetto (niente
rotazione terrestre). Esplorate e scartate con evidenza numerica due
alternative (obiettivo auto-consistente, tempo di bruciamento libero)
prima di arrivare a questa spiegazione — nessuna scorciatoia, il
percorso è documentato per intero in questo log. Prossimo step
proposto: Step 10 (validazione e documentazione dei limiti — questa
scoperta e la sua chiusura sono un candidato naturale per VALIDATION.md).

---

### 2026-08-18 — Ciclo 11 (io, pipeline snellita): piano Step 10 — VALIDATION.md

**Step:** 10 — Validazione e documentazione dei limiti (VALIDATION.md,
confronto concettuale con i progetti di riferimento)
**Stato:** pianificazione completata (io, senza sub-agente planner, per
lo stesso motivo già seguito dagli step precedenti: step di
documentazione, nessuna nuova fisica).

#### 1. Contesto

Step 1-9 + Ciclo 10 completi, 91/91 test verdi, commit locale `e1cb826`.
Questo step non introduce nuova fisica: raccoglie in un documento unico
(`VALIDATION.md`, nuovo file in root) tutte le approssimazioni
dichiarate in CLAUDE.md e i loro impatti REALI già scoperti e
quantificati durante il progetto (nessuna nuova stima), più un
confronto concettuale onesto con i tre progetti di riferimento citati
in CLAUDE.md.

#### 2. Ricerca preliminare (WebFetch, letto ma non copiato)

Ricerca fattuale sui tre repository citati in CLAUDE.md come
"riferimento concettuale (NON copiare codice)":

- **axelstr/gravity_turn_simulation**: 2D, MULTISTADIO, guida attiva
  reale (cutoff automatico all'apogeo target + manovra di
  circolarizzazione) — più empirica/euristica della tangente lineare
  di questo progetto (derivata analiticamente da Perkins/PEG). Nessuna
  validazione con dati reali dichiarata nel repo.
- **bvermeulen/Rocket-and-gravity-turn**: 2D, MONOSTADIO, ha ANCHE
  controllo ottimo della spinta verso obiettivi orbitali (risolutore
  CasADi, ottimizzazione numerica diretta/collocazione) — sofisticazione
  comparabile alla guida di questo progetto, ma con metodo diverso
  (collocazione numerica vs shooting analitico). Nessuna validazione
  con dati reali dichiarata.
- **b-adkins/pyrocket**: attualmente 1D (2D è roadmap, non implementato),
  MULTISTADIO, nessuna guida attiva sofisticata (gravity turn "in
  roadmap"). Orientato a scopo didattico/stile KSP.

Questo AGGIORNA (rende più preciso, non contraddice) il claim originale
di CLAUDE.md "i progetti individuali reali... usano tutti punto
materiale 2D": resta vero per la dinamica, ma va precisato che due dei
tre hanno comunque guida attiva sofisticata (uno multistadio con
cutoff automatico, l'altro con controllo ottimo via CasADi) — nessuno
dei tre usa pero' una legge di guida analiticamente derivata e citata
in letteratura (Perkins/PEG) come questo progetto, ne' dichiara un
confronto quantitativo con un benchmark reale con margine esplicito.

#### 3. Design — struttura di VALIDATION.md

Nuovo file `VALIDATION.md` in root, 5 sezioni:

1. **Sintesi del check di validazione obbligatorio**: Δv ideale totale
   9138.15 m/s, dentro il range 9.1-10.0 km/s, margine 38 m/s sopra il
   limite inferiore — riassunto, rimanda a STATUS.md Ciclo 7/10 per il
   dettaglio del calcolo (non ripetuto per esteso).
2. **Tabella delle approssimazioni dichiarate in CLAUDE.md, con impatto
   MISURATO (non stimato ora) dove disponibile:**
   - *Atmosfera esponenziale*: errore ~26% a 11km, ~2.5% a 25km rispetto
     al modello di riferimento (Step 1, test_atmosfera.py).
   - *Cd costante*: nessun impatto quantificato disponibile nel
     progetto (nessun dato di confronto Cd-vs-Mach) — dichiarato
     onestamente come limite non misurato, non presentato come
     "presumibilmente trascurabile".
   - *Gravity turn senza sollievo centripeto*: 12.2%-17.3% di g alla
     velocità di fine Stadio 1 (Ciclo 7 Fase A) — NON corretto (solo la
     fase esoatmosferica lo è, Ciclo 7 Fase A), asimmetria accettata e
     dichiarata esplicitamente, non risolta.
   - *Niente rotazione terrestre*: bonus escluso ≈408.74 m/s
     (`465.1*cos(28.5°)`, sito equatoriale/subtropicale) — **dimostrato
     essere la spiegazione quantitativa completa del deficit di perigeo
     scoperto al Ciclo 10** (deficit ~67 m/s, coperto dal bonus con
     margine). Raccontato per esteso qui come l'esempio più forte del
     progetto di "approssimazione dichiarata → impatto misurato →
     conseguenza verificata due volte in modo indipendente", non solo
     linkato a STATUS.md.
   - *Drag inferiore al range di letteratura*: 16.7 m/s di perdita
     ideale da drag, contro 100-400 m/s tipicamente citati in
     letteratura generica — spiegato (non un errore nascosto) via
     confronto Max-Q con dati reali Falcon 9 (Ciclo 9/10,
     critic-ingegnere): il profilo di traiettoria di questo progetto
     attraversa la zona di massima pressione dinamica più rapidamente/
     con Cd*A minore del velivolo reale completo.
3. **Fenomeno delle radici multiple nel problema inverso tangente
   lineare** (Step 6/8): non un'approssimazione fisica ma una
   caratteristica nota del metodo numerico (shooting via
   scipy.optimize) — documentato come limite metodologico: richiede un
   guess nel bacino di attrazione corretto, nessuna garanzia di
   convergenza globale. Mitigato (non eliminato) dal safeguard a doppio
   guess costruito allo Step 6 e riusato allo Step 8.
4. **Confronto concettuale con i progetti di riferimento** (tabella,
   dati dalla ricerca del punto 2, parafrasati non citati/copiati):
   scope (2D/3D), multistadio, tipo di guida, validazione con dati
   reali dichiarata. Posiziona onestamente questo progetto (guida
   analiticamente derivata e citata in letteratura, validazione
   quantitativa con margine esplicito e limiti dichiarati) senza
   sminuire gli altri (approcci diversi, non necessariamente inferiori
   — es. il controllo ottimo CasADi di bvermeulen è un metodo valido,
   solo diverso da quello scelto qui).
5. **Estensioni future dichiarate**: riprese da CLAUDE.md, con lo stato
   di avanzamento reale — es. lo Step 8 (staging in tangente lineare)
   copre già qualcosa che la formulazione originale del vincolo
   "Multistadio" non specificava esplicitamente per la fase
   esoatmosferica; le altre (rotazione terrestre, 6-DOF, Cd variabile,
   ottimizzazione della traiettoria) restano fuori scope come
   dichiarato, senza modifiche.

#### 4. Esecuzione (pipeline snellita)

1. Questo piano, scritto qui in STATUS.md.
2. **reviewer**: verifica accuratezza (nessun numero sovra/sottostimato
   rispetto a STATUS.md) e che il confronto con i progetti di
   riferimento sia onesto e non copi testo/codice (solo parafrasi
   fattuale).
3. Io scrivo `VALIDATION.md` (documentazione pura, nessun `coder`
   necessario).
4. **critic-ingegnere**: verifica che ogni numero citato in
   VALIDATION.md corrisponda esattamente a un valore già verificato in
   STATUS.md (nessuna nuova stima non tracciata).
5. Pulizia, aggiornamento STATUS.md, report — fatti da me.
6. Mi fermo e attendo conferma dell'utente prima di passare allo
   Step 11 (pulizia, README, preparazione per GitHub).

#### 5. File coinvolti

- `VALIDATION.md` (nuovo, root del progetto)
- `STATUS.md` (questo piano)
- Nessuna modifica al codice Python, nessun nuovo test (step di sola
  documentazione)

**Prossimo:** dispatch del sub-agente reviewer su questo piano.

#### 6. Esito ciclo (reviewer + verifica indipendente + critic-ingegnere)

- **reviewer**: cross-check numerico di tutti i valori del piano contro
  STATUS.md — nessuna discrepanza trovata. Segnalati 6 punti di
  miglioramento, di cui 2 vincolanti (requisiti già imposti in
  precedenza nel progetto, non solo suggerimenti nuovi): (1) mancava il
  numero di sensibilità Isp-SL (8763.81 m/s) nella sintesi del check
  Δv — obbligo esplicito del reviewer al Ciclo 7 ("riportare
  ESPLICITAMENTE ANCHE il numero con Isp SL... come limite inferiore
  di sensibilità"); (2) nessuna verifica indipendente delle
  affermazioni sui tre progetti esterni (un'unica ricerca WebFetch
  fatta da un solo agente, incoerente con lo standard di doppia
  verifica già applicato ad ogni scoperta importante del progetto).
  Altri 4 suggerimenti minori: citazione URL/data per la ricerca sui
  repo esterni, raccordo esplicito Δv-ideale-vs-orbita-raggiunta,
  asimmetria non dichiarata guida-tangente-lineare-derivata-per-g-costante
  vs dinamica-integrata-con-g-variabile, numeri di supporto per la
  Sezione 3 (radici multiple).
- **Verifica indipendente aggiuntiva (io)**: dispatch di un secondo
  agente (general-purpose, con WebFetch) per rifare la ricerca sui tre
  repository esterni SENZA fornirgli le mie conclusioni precedenti,
  per un controllo davvero indipendente (risponde al punto 2 del
  reviewer). Risultato: conferma sostanziale di tutti i fatti già
  raccolti, con una precisazione (axelstr non dichiara esplicitamente
  la dimensionalità 2D nel README — è un'inferenza dal contesto, non
  un fatto dichiarato — corretto in VALIDATION.md per non
  sovra-affermare).
- **Scrittura VALIDATION.md (io)**: tutti i 6 punti del reviewer
  incorporati: aggiunta la sensibilità Isp-SL in Sezione 1 con nota
  esplicativa sul perché Isp-vuoto è la scelta di base ma il margine è
  sensibile a quella scelta; aggiunto un paragrafo di raccordo che
  distingue esplicitamente "capacità ideale del veicolo" (check Δv)
  da "stato cinematico finale raggiunto" (orbita); aggiunta una nuova
  sottosezione sull'asimmetria guida-tangente-lineare/dinamica
  integrata; aggiunto l'esempio numerico concreto del Ciclo 6
  (A=0.258,B=+0.002 vs A_vero=0.3,B_vero=-0.002) alla Sezione 3;
  aggiunti i link diretti ai tre repository nella tabella di
  confronto Sezione 4, con nota sulla dimensionalità di axelstr
  corretta a "non dichiarata esplicitamente, verosimilmente...".
- **critic-ingegnere**: verifica puntuale di tutti i 13 numeri citati
  in VALIDATION.md contro STATUS.md — tutti CONFERMATI (con
  arrotondamenti dichiarati e coerenti, es. 16.73→16.7 m/s,
  -22.77→-22.8 km), nessun numero nuovo non tracciato o alterato.
  Confermata anche la coerenza della Sezione 4 (confronto progetti
  esterni) con il piano di questo ciclo. **Verdetto: il documento è
  pronto così com'è.**

**Ciclo 11: COMPLETATO.** `VALIDATION.md` creato in root, 5 sezioni
come da piano, tutti i numeri tracciabili a STATUS.md, nessuna nuova
fisica introdotta, confronto con i progetti di riferimento basato su
ricerca reale e verificato due volte in modo indipendente. Nessuna
modifica al codice Python, nessun nuovo test (step di sola
documentazione, coerente con la sua natura). Prossimo step proposto:
Step 11 (pulizia, documentazione, README, preparazione per GitHub).
