# Simulazione di ascesa e guida di un lanciatore multistadio

Simulazione a punto materiale (2D, piano verticale) dell'ascesa
propulsa di un lanciatore multistadio, con **guida attiva** (non
traiettoria scriptata) verso un'orbita target: gravity turn in fase
atmosferica, guida a tangente lineare (linear-tangent steering,
Perkins/PEG) in fase esoatmosferica. Validata contro dati pubblici di
un lanciatore reale (Falcon 9).

Lo scope (2D, punto materiale, senza rotazione terrestre) è una scelta
dichiarata e motivata — vedi [CLAUDE.md](CLAUDE.md) per il confronto
con i progetti open source di riferimento e la tabella completa dei
vincoli tecnici.

## Caratteristiche principali

- **Gravity turn** (fase atmosferica): allineamento spinta-velocità
  dopo un kick angle iniziale (Culler & Fried, 1957).
- **Guida a tangente lineare** (fase esoatmosferica): `tan(θ) = A + B·t`,
  ottimo matematico per la sola spinta secondo il principio del
  massimo di Pontryagin (Perkins) — la stessa base teorica del PEG
  (Powered Explicit Guidance) usato dallo Space Shuttle.
- **Multistadio**, con eventi di staging (cambio massa discontinuo)
  gestiti come eventi terminali dell'integrazione ODE — anche a metà
  della fase di guida a tangente lineare, con ricalcolo dei
  coefficienti A/B al momento dello staging.
- **Problema inverso**: root-finding di A/B per un target di velocità
  terminale, con un safeguard esplicito contro il fenomeno delle
  radici multiple del sistema (vedi [VALIDATION.md](VALIDATION.md)).
- **Validazione quantitativa** contro specifiche pubbliche Falcon 9
  Block 5, con verifica di meccanica orbitale (perigeo/apogeo) sullo
  stato finale raggiunto, non solo un check di delta-v aggregato.
- **Visualizzazione**: grafici quota/velocità nel tempo, traiettoria,
  animazione della traiettoria in due versioni — una semplice (punto +
  linea) e una stilizzata, con un simbolo del razzo orientato in tempo
  reale lungo la direzione di volo reale (ruota visibilmente durante il
  gravity turn), scia con dissolvenza colorata per stadio, e un cambio
  di colore/dimensione esatto all'istante di stage separation. Nessuna
  fisica nuova in nessuna delle due: puro rendering di numeri già
  calcolati e validati.

## Struttura del repository

```
lanciatore/
  costanti.py               costanti fisiche (G0, R_TERRA, MU_TERRA, RHO0, H_SCALA, CD)
  atmosfera.py               modello di densità atmosferica esponenziale
  gravita.py                 accelerazione di gravità in funzione della quota
  dinamica.py                dinamica punto materiale, ascesa verticale pura (caso di test)
  guida.py                   gravity turn (fase atmosferica)
  guida_esoatmosferica.py    guida a tangente lineare + risoluzione del problema inverso A/B
  staging.py                 multistadio in fase di gravity turn
  staging_esoatmosferico.py  staging durante la guida a tangente lineare
  validazione.py             caso di validazione Falcon 9, delta-v, meccanica orbitale
  visualizzazione.py         grafici e animazione della traiettoria

tests/                       100 test, uno o più per modulo
STATUS.md                    log dettagliato di sviluppo, ciclo per ciclo
VALIDATION.md                limiti del modello, impatti misurati, confronto con progetti di riferimento
CLAUDE.md                    scope del progetto, vincoli tecnici, regole di lavoro
```

## Installazione

```bash
python -m venv venv
venv\Scripts\activate   # su Windows; su Unix: source venv/bin/activate
pip install -r requirements.txt
```

## Eseguire i test

```bash
pytest
```

100 test coprono ogni modulo, incluse verifiche di regressione contro
valori noti (formule chiuse, oracoli di meccanica orbitale) e casi
espliciti che documentano i limiti noti come "spiegati", non aperti.

## Generare la validazione e i grafici

```bash
python -m lanciatore.visualizzazione
```

Rigenera in `output/` (non tracciato in git) i grafici quota/velocità
nel tempo, la traiettoria, l'animazione semplice (`traiettoria.gif`) e
quella stilizzata (`traiettoria_avanzata.gif`, razzo orientato + scia +
stage separation), a partire dal risultato di
`lanciatore.validazione.assembla_traiettoria()` — nessun numero è
hardcodato nei grafici. L'angolo di orientamento del razzo nella
versione avanzata (`lanciatore.visualizzazione.angolo_volo_gradi`) è
anch'esso derivato da stato già integrato (gamma del gravity turn,
`atan2(vh, vx)` della tangente lineare), non una nuova grandezza fisica.

## Risultato principale

Il caso di validazione (parametri Falcon 9 Block 5, dati pubblici) dà
un **delta-v ideale totale di 9138.15 m/s**, dentro il range di
plausibilità 9.1-10.0 km/s per un'orbita LEO, con un margine di 38 m/s
(0.42%) sopra il limite inferiore. Dettaglio completo, inclusa la
sensibilità a diverse scelte di Isp e il confronto meccanica-orbitale
sullo stato finale raggiunto (non solo il delta-v aggregato), in
[VALIDATION.md](VALIDATION.md) e [STATUS.md](STATUS.md).

## Limiti noti

Ogni approssimazione dichiarata nello scope del progetto (atmosfera
esponenziale, Cd costante, niente sollievo centripeto nel gravity
turn, niente rotazione terrestre) ha un impatto misurato e discusso
per esteso in **[VALIDATION.md](VALIDATION.md)**, insieme al fenomeno
delle radici multiple nel problema inverso e a un confronto onesto con
tre progetti open source di riferimento concettuale. Non ripetuto qui.

## Documentazione

- [CLAUDE.md](CLAUDE.md) — scope del progetto, vincoli tecnici non
  negoziabili, regole di lavoro.
- [STATUS.md](STATUS.md) — log completo di sviluppo, ciclo per ciclo,
  incluse le decisioni prese e il perché.
- [VALIDATION.md](VALIDATION.md) — limiti del modello, impatti
  misurati, confronto con progetti di riferimento.

## Licenza

[MIT](LICENSE)
