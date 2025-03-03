import argparse
import timeit

from crispr import crisprcas, crisprcas_next_to

startTime = timeit.default_timer()

parser = argparse.ArgumentParser(description="Diseño de guias CRISPR-Cas9 o Cas12 a partir de una secuencia en formato FASTA",
                                 epilog="Ejemplos de uso:\n "
                                        "python run_crispr.py --function crisprcas --fasta ../data/editar.fasta --cas 9 --section all --graph"
                                        "python run_crispr.py --function crisprcas_next_to --fasta ../data/editar.fasta --cas 12 --section_nt 500"
                                        "NOTA: Los archivos se guardan en la carpeta 'results'")
parser.add_argument("--function",
                        choices=["crisprcas", "crisprcas_next_to"],
                        default="crisprcas",
                        help="Tipo de guias CRISPR deseadas: crisprcas, para secciones en general, o crisprcas_next_to, para guias cercanas a un nucleotido específico (valor por defecto: crisprcas)")
parser.add_argument("--fasta",
                        required=True,
                        help="Ruta al archivo FASTA que contiene el gen o la secuencia de nucleotidos a editar")
parser.add_argument("--cas",
                        type=int,
                        choices=[9, 12],
                        default=9,
                        help="Tipo de sistema CRISPR deseado: 9 o 12 (valor por defecto: 9)")
parser.add_argument("--section",
                        choices=["all", "start", "middle", "end"],
                        default="all",
                        help="Seccion del gen en la que se desean identificar las guias: all, start, middle o end (valor por defecto: all)")
parser.add_argument("--section_nt",
                        type=int,
                        default=0,
                        help="Posicion del gen en la que se desean identificar las guias ±100 nucleotidos (valor por defecto: 0, es decir, las guias dentro de los primeros 100 nucleotidos)")
parser.add_argument("--graph",
                        type=bool,
                        default=False,
                        help="Indicar si se desea generar graficos de las guias: True o False (valor por defecto: False)")


args = parser.parse_args()

if args.function == "crisprcas":
    result = crisprcas(args.fasta, args.cas, args.section, args.graph)
else:
    result = crisprcas_next_to(args.fasta, args.cas, args.section_nt, args.graph)

endTime = timeit.default_timer()

print(f"Fin del proceso "
      f"Tiempo de ejecución {endTime - startTime} segundos")