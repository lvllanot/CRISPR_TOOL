# CRISPR Guide Designer

Este proyecto es una herramienta en Python que permite diseñar y seleccionar las mejores guías CRISPR a partir de un archivo FASTA con una secuencia de nucleótidos. El programa procesa la secuencia, genera candidatos de guías (para sistemas CRISPR Cas9 o Cas12), verifica su especificidad y genera una base de datos (TSV) con las guías seleccionadas. Además, puedes ejecutar el programa desde la línea de comandos y elegir si deseas generar gráficos que ilustren la posición de las guías a lo largo de la secuencia.

---
## Procesamiento

- **Diseño de candidatos de guías:**  
  Se diseñan candidatos de guías de 20 nucleótidos según el sistema CRISPR especificado (Cas9 o Cas12).

- **Verificación de especificidad:**  
  Se utiliza BLAST local para verificar la especificidad de las guías.

- **Filtrado de guías:**  
  Se filtran las guías basándose en criterios como el contenido de GC y el número de alineamientos.

- **Salida:**  
  Se genera una base de datos en formato TSV con las mejores guías y sus características (orientación, posición, contenido de GC, etc.).

- **Opciones de visualización:**  
  Posibilidad de generar gráficos que muestren la distribución de las guías sobre la secuencia.

- **Ejecución desde línea de comandos:**  
  Utiliza `argparse` para definir y gestionar los argumentos de entrada.

---
## Requisitos

- **Python 3.x**
- **Pandas**
- **Matplotlib**
- **BLAST+** (instalado y en el PATH)
- **Seaborn**

---
## Ejemplo

**Ejemplo básico**
Para ejecutar el programa usando la función crisprcas (guías para secciones generales):

```bash
python run_crispr.py --function crisprcas --fasta "../data/editar.fasta" --cas 9 --section all --graph True

```
** Ejemplo para guías cerca de un nucleótido específico **
Para ejecutar la función crisprcas_next_to y especificar la posición (número entero) donde se deben buscar guías:

```bash
python run_crispr.py --function crisprcas_next_to --fasta ../data/editar.fasta" --cas 9 --section_nt 500 --graph False
```

---
## Argumentos

**--function:**

- **Opciones:** crisprcas, crisprcas_next_to  
- **Descripción:** Define el modo de funcionamiento (por defecto crisprcas).  
  - **crisprcas:** Diseña guías para secciones generales del gen.  
  - **crisprcas_next_to:** Diseña guías cercanas a un nucleótido específico.

**--fasta:**

- **Descripción:** Ruta al archivo FASTA del gen o secuencia a editar.

**--cas:**

- **Opciones:** 9 o 12  
- **Descripción:** Define el tipo de sistema CRISPR (por defecto 9).

**--section:**

- **Opciones:** all, start, middle, end  
- **Descripción:** Seccion del gen en la que se desean identificar las guias (solo para crisprcas).

**--section_nt:**

- **Tipo:** entero  
- **Descripción:** Posición del gen para considerar guías cercanas (solo para crisprcas_next_to).

**--graph:**

- **Descripción:** Indicar si se desea generar graficos de las guías(por defecto False).

---
## Resultados

Los archivos generados (por ejemplo, la base de datos TSV y los gráficos) se guardan en la carpeta `results`.

---
## Notas

- Asegúrate de que BLAST+ esté instalado y configurado en tu sistema para la verificación de las guías.

---
## Licencia

Este proyecto está licenciado bajo la [MIT License](LICENSE).
