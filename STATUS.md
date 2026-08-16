# Stato del progetto — aggiornato ad ogni ciclo

## Step
- [x] Step 1: Setup ambiente, modello atmosferico esponenziale, costanti fisiche
- [x] Step 2: Dinamica a punto materiale (spinta, gravita', drag) — singolo stadio, ascesa verticale pura, senza guida (caso di test piu' semplice)
- [x] Step 3: Gravity turn (fase atmosferica)
- [ ] Step 4: Guida a tangente lineare (fase esoatmosferica) — verifica contro derivazione analitica nota
- [ ] Step 5: Multistadio (eventi di staging, cambio massa discontinuo)
- [ ] Step 6: Caso di validazione con dati reali di un lanciatore pubblico + confronto delta-v vs benchmark ~9.1-10.0 km/s
- [ ] Step 7: Visualizzazione (traiettoria numerica + animazione)
- [ ] Step 8: Validazione e documentazione dei limiti (VALIDATION.md, confronto concettuale con i progetti di riferimento)
- [ ] Step 9: Pulizia, documentazione, README, preparazione per GitHub

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
