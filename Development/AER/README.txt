Le funzioni sono in ordine di accuratezza:
1) class1 --> accuratezza bassa e geometria semplice
2) class2 --> accuratezza media e geometria piu elaborata (possibilità di inserire ali e coda)
3) class3 --> accuratezza alta, ottimizzazione aerodinamica sarà basatasu questa funzione


INPUT FUNZIONI:
- Geometria (definita all'interno della funzione)
- h: altitudine
- Mach
- Angolo d'attacco

OUTPUT FUNZIONI:
- Cd, Cl (Ca, Cn)
- Cm e Xcp (solo class2 e class3)

OSS: le funzioni ready to run sono quelle: " _fun "
