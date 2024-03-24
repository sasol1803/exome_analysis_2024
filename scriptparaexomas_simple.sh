#Pipeline de análisis de exomas
#Ejemplo para análisis de un caso individual, cambiar nombre de archivo inicial según la muestra a estudiar. 

#Antes de empezar creamos un enviroment de conda con los programas que necesitamos
conda create -n exome python=3.8
conda activate exome

#Instalación de softwares en enviroment de trabajo

conda install -c bioconda fastqc
conda install -c bioconda fastp
conda install -c bioconda multiqc
conda install -c bioconda bwa 
conda install -c bioconda samtools 
conda install -c bioconda qualimap
conda install -c bioconda picard
conda install -c bioconda freebayes
conda install -c bioconda vcftools
conda install -c bioconda bcftools
conda install -c bioconda igv

#1. Análisis primario de calidad de las secuencias.
#1.1 Análisis de calidad de los archivos crudos con fastqc

#Partimos del directorio donde están los fastq.gz de nuestro exoma. Tenemos R1 y R2 por muestra.
cd home/laura/exomas/exome_01

#Hacemos el qc con fastqc y movemos los archivos de salida a una nueva carpeta
fastqc *.fastq.gz
mkdir FastQC_raw_reads
mv *fastqc* FastQC_raw_reads

#Nos vamos a la carpeta con el output y vemos los .html de salida
cd FastQC_raw_reads
firefox *.html

#1.2 Limpieza de las secuencias: eliminar adaptadores y trimming

#Volvemos a la carpeta exome_01 y ejecutamos fastp

cd ..
fastp -i indice_R1_001.fastq.gz -I indice_R2_001.fastq.gz -o indice1.clean.fastq.gz -O indice2.clean.fasq.gz --trim_poly_g --detect_adapter_for_pe --cut_front 25 --cut_tail 25 --cut_mean_quality 25 -l 60 -h report_fastp.html

#Repetimos el control de calidad
fastqc *clean.fastq.gz

#Creamos una nueva carpeta para guardar los resultados y observamos de nuevo los .html
mkdir FastP_result
mv *clean* *fastp* FastP_result/
cd FastP_result
firefox *fastqc*

#Opcional, ver calidad con multiqc, en la carpeta de los archivos limpios
multiqc .

#Creamos una carpeta con los fastq iniciales y otra con los fastq clean
cd..
mkdir fastq_raw 
mkdir fastq_clean
mv *fastq.gz fastq_raw
mv FastP_result/*clean* fastq_clean

#2. Mapeo/Alineamiento a genoma de referencia
#2.1 Localizar el genoma de referencia
#Indexar el genoma (solamente es necesario en su primer uso)

hg38='/home/biodatabases/GRCh38.p14/GRCh38.p14.fna'
bwa index $hg38

#Conversión de archivos limpios a .sam
mkdir Mapeo
cd Mapeo
bwa mem -a $hg38 ../fastq_clean/indice1.clean.fastq.gz ../fastq_clean/indice2.clean.fastq.gz -o indice.sam 2> indice.out

#2.2 Conversión de archivos .sam a .bam y procesamiento
samtools view -bS indice.sam > indice.bam

#Visualización y visualización por calidad
samtools view -h indice.bam | more
samtools view -q 30 indice.bam | more

#Ordenar e indexar
samtools sort indice.bam > indice.sorted.bam
samtools index indice.sorted.bam

#2.3 Visualización del archivo .bam con IGV 

#2.4 Análisis de la calidad del BAM

#Visualización de estadísticas básicas
samtools flagstat indice.sorted.bam > indice.flagstat.txt

#creación de informe da calidad con quialimap
qualimap bamqc -bam indice.sorted.bam -c -nt 2 -gd HUMAN -gff genome.gff -outformat html -java-mem-size=32G

#3. Identificación de variantes
#3.1 Identificación de duplicados con PicardTools

picard MarkDuplicates --INPUT indice.sorted.bam --OUTPUT indice.dedup.sorted.bam --METRICS_FILE markDuplicatesMetrics.txt --ASSUME_SORTED True

#Añadimos readgroups 
picard AddOrReplaceReadGroups -I indice.dedup.sorted.bam -O indice.dedup.sorted.RG.bam -RGID M02899 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM S8

#Volvemos a indexar el bam
samtools index indice.dedup.RG.sorted.bam

#Podemos visualizar el archivo en IGV

#3.2 Identificación de variantes
#Con freebayes
#25 Minimo de lecturas a mapear en una posición para considerarse alelo ALT #Poner genoma 38 completo 

freebayes -C 25 -f $hg38 indice.dedup.RG.sorted.bam > indice.vcf

#Estadísticas
rtg vcfstats indice.vcf > indice.vcfstats

#3.3 Filtrado de variantes

#SNPs
vcftools --vcf indice.vcf --keep-only-indels --recode --recode-INFO-all --out indice_indels.vcf
#INDELS
vcftools --vcf indice.vcf --remove-indels --recode --recode-INFO-all --out indice_snvs.vcf


#Alternativa con GATK
gatk='/home/bioinfotools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar'
threads='8'
RAM='Xmx32g'
illumina_exome='/home/biodatabases/bedpanels/gencode40_refseq_downpadding30.bed'

java -$RAM -jar $gatk HaplotypeCaller --native-pair-hmm-threads $threads -R $hg38 -I indice.dedup.RG.sorted.bam -O indice.vcf.gz -L $illumina_exome
java -$RAM -jar $gatk SortVcf \
	-I indice.vcf.gz \
	-O indice.sorted.vcf.gz

# SNP
java -jar $gatk SelectVariants -R $hg38 -V indice.sorted.vcf.gz -select-type SNP -O indice.sorted.snp.vcf.gz
java -$RAM -jar $gatk VariantFiltration -R $hg38 -V indice.sorted.vcf.gz \
-filter "QUAL < 30" \
--filter-name "Quality" \
-filter "DP < 20" \
--filter-name "LowDepth" \
-O indice.sorted.snp.filtered.vcf.gz

# INDEL
java -jar $gatk SelectVariants -R $hg38 -V indice.sorted.vcf.gz -select-type INDEL -O indice.sorted.indels.vcf.gz
java -jar $gatk VariantFiltration -R $hg38 -V indice.sorted.vcf.gz \
-filter "QUAL < 30" \
--filter-name "Quality" \
-O indice.sorted.indels.filtered.vcf.gz


#3.4 Visualización de las variantes con IGV #Cargamos sondas, BAM y VCF

#4. Anotación de variantes
#Partimos del .vcf con todas las variantes

#Anotación con ANOVAR
annovar='/home/bioinfotools/annovar'
chrnames='/home/biodatabases/bedpanels/refseq_modified.dict'

bcftools view -r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,\
NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,\
NC_000023.11,NC_000024.10,NC_012920.1 -o indice.sorted.vcf.gz -O z indice.sorted.snp.filtered.vcf.gz && tabix indice.clean.vcf.gz

bcftools annotate --rename-chr $chrnames -o indice.chr.vcf.gz -O z indice.clean.vcf.gz && tabix indice.chr.vcf.gz

perl $annovar/table_annovar.pl indice.chr.vcf.gz $annovar/humandb/ --buildver hg38 \
--out --remove --protocol knownGene,refGeneWithVer,avsnp150,clinvar_20220320,cosmic70,esp6500siv2_all,1000gALL,1000gEUR,exac03,gnomad_exome,ljb26_all \
--operation g,g,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --polish
bgzip indice.hg38.anno.vcf

#Alternativa de anotación vía web con wANNOVAR o VEP. 

#Organización final de directorios
#Dentro de la carpeta de Mapeo creamos carpetas con los diferentes archivos que hemos creado en pasos anteriores
mkdir SAM
mkdir BAM
mkdir flagstat
mdir Qualimap
mkdir Markduplicated
mkdir VCF

#movemos los archivos finales a sus carpetas y eliminamos archivos intermediarios no necesarios

