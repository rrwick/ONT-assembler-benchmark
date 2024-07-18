This file contains my commands for selecting and preparing genomes from the two big P2 runs of April 2024.

I did most of this analysis on Roosta, except the basecalling was done on Spartan via SLURM.



# Download and prep Illumina reads

Tools:
* fastp v0.23.4

Download reads:
```bash
cd ~/2024-06_new_assembler_benchmark
mkdir -p illumina_reads; cd illumina_reads
curl https://bioinformatics.mdu.unimelb.edu.au/~cwwalsh/tmp/D02AIQSG2JDRI7Y12OK8VJ0VD/download.sh | bash
curl https://bioinformatics.mdu.unimelb.edu.au/~cwwalsh/tmp/ZT0DP8GLGHW7CSPYXHXGBUOOC/download.sh | bash
curl https://bioinformatics.mdu.unimelb.edu.au/~bisognor/tmp/2M9WUWK282FTMA8I0UPUW4FP6/download.sh | bash
curl https://bioinformatics.mdu.unimelb.edu.au/~bisognor/tmp/6RQLWJURHQKV2H225AISVOMO2/download.sh | bash

cd ~/2024-06_new_assembler_benchmark
for id in $(ls illumina_reads/*_R1*.fastq.gz | grep -oP "AUSMDU\d+"); do
    mkdir -p "$id"/reads
    mkdir -p "$id"/reads_qc
    mv illumina_reads/"$id"_*R1*.fastq.gz "$id"/reads/illumina_1.fastq.gz
    mv illumina_reads/"$id"_*R2*.fastq.gz "$id"/reads/illumina_2.fastq.gz
done
rmdir illumina_reads
```

Get some read stats and check that paired files contain the same read count:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"/reads
    fast_count illumina_1.fastq.gz | awk '{print $2 " " $3}' | read count1 bases1
    fast_count illumina_2.fastq.gz | awk '{print $2 " " $3}' | read count2 bases2
    if [[ "$count1" -ne "$count2" ]]; then
        echo "Error: Count mismatch in $id. count1: $count1, count2: $count2"
        break
    fi
    printf "$id\t$count1\t$bases1\t$bases2\n"
done
```

QC:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"/reads_qc
    fastp --in1 ../reads/illumina_1.fastq.gz --in2 ../reads/illumina_2.fastq.gz --out1 illumina_1.fastq.gz --out2 illumina_2.fastq.gz --unpaired1 illumina_u.fastq.gz --unpaired2 illumina_u.fastq.gz
done
```


And get post-QC read stats:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"/reads_qc
    fast_count illumina_1.fastq.gz | awk '{print $2 " " $3}' | read count1 bases1
    fast_count illumina_2.fastq.gz | awk '{print $2 " " $3}' | read count2 bases2
    if [[ "$count1" -ne "$count2" ]]; then
        echo "Error: Count mismatch in $id. count1: $count1, count2: $count2"
        break
    fi
    printf "$count1\t$bases1\t$bases2\n"
done
```





# Illumina-only assembly

Even though I'm not going to be using Illumina-only assemblies, they are useful for QCing the Illumina read sets, i.e. helping me to exclude samples with bad Illumina reads.

Tools:
* Unicycler v0.5.1
* SPAdes v4.0.0
* GFA-dead-end-counter v0.1.0

Assemble with Unicycler:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    if [ -f unicycler_illumina/assembly.gfa ]; then continue; fi
    unicycler -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -o unicycler_illumina -t 48
done
```

Get some assembly stats:
```bash
cd ~/2024-06_new_assembler_benchmark
for gfa in AUSMDU*/unicycler_illumina/assembly.gfa; do
    id=$(echo "$gfa" | grep -oP "AUSMDU\d+")
    fasta=${gfa/gfa/fasta}
    depth_filter_gfa=${gfa/assembly/002_depth_filter}

    best_k=$(cat "$depth_filter_gfa" | grep -P "^L" | grep -oP "\d+M" | uniq | sed 's/M//' | xargs printf "%03d\n")
    best_k_gfa=${gfa/assembly/001_spades_graph_k"$best_k"}

    best_k_gfa_size=$(any2fasta "$best_k_gfa" 2> /dev/null | seqkit stats -T | tail -n1 | cut -f5)
    depth_filter_gfa_size=$(any2fasta "$depth_filter_gfa" 2> /dev/null | seqkit stats -T | tail -n1 | cut -f5)
    contamination_size=$[$best_k_gfa_size - $depth_filter_gfa_size]

    final_contigs=$(fast_count "$fasta" | cut -f2)
    final_size=$(fast_count "$fasta" | cut -f3)
    estimated_genome_size=$(scripts/estimate_genome_size.py "$gfa")
    deadends=$(deadends "$gfa")

    printf "$id\t$final_contigs\t$final_size\t$estimated_genome_size\t$contamination_size\t$deadends\n"
done
```
The estimated genome size (via my script) could sometimes be an overestimate (due to tangled small plasmids), but it is still better than the estimate I'd get from Bandage.




# Prepare pod5 files

Make directories (on Spartan):
```bash
cd /data/scratch/projects/punim1894
mkdir O2024-029 O2024-030
mkdir O2024-030/pod5
```

Transfer reads (on Roosta):
```bash
cd /home/damg/data/O2024-029
scp -r pod5 spartan:/data/scratch/projects/punim1894/O2024-029
cd /home/damg/data/O2024-030
scp *.pod5 spartan:/data/scratch/projects/punim1894/O2024-030/pod5
```

Count reads (on Spartan):
```bash
cd /data/scratch/projects/punim1894/O2024-029/pod5
pod5 view *.pod5 | tail -n+2 | wc -l
cd /data/scratch/projects/punim1894/O2024-030/pod5
pod5 view *.pod5 | tail -n+2 | wc -l
```

Pod5 read counts:
* O2024-029: 23767065
* O2024-030: 15823549



# Simplex basecalling

Tools:
* Dorado v0.7.0

Dorado commands (on Spartan):
```bash
export PATH=/home/rrwick/programs/dorado-0.7.0-linux-x64/bin:"$PATH"

cd /data/scratch/projects/punim1894/O2024-029
dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam
# sbatch --job-name=O2024-029 --time=100:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-a100 --gres=gpu:1 --wrap "~/programs/dorado-0.7.0-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam"
dorado summary reads.bam > summary.tsv
dorado demux --output-dir reads --no-classify reads.bam
cd reads
for b in *.bam; do
    samtools fastq -T '*' "$b" | tr '\t' ' ' | paste - - - - | sort | tr '\t' '\n' > ${b/bam/fastq}
done
for f in SQK*.bam SQK*.fastq; do
    new_f=${f/SQK-RBK114-96_/}
    mv "$f" "$new_f"
done

cd /data/scratch/projects/punim1894/O2024-030
dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam
# sbatch --job-name=O2024-030 --time=70:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-a100 --gres=gpu:1 --wrap "~/programs/dorado-0.7.0-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam"
dorado summary reads.bam > summary.tsv
dorado demux --output-dir reads --no-classify reads.bam
cd reads
for b in *.bam; do
    samtools fastq -T '*' "$b" | tr '\t' ' ' | paste - - - - | sort | tr '\t' '\n' > ${b/bam/fastq}
done
for f in SQK*.bam SQK*.fastq; do
    new_f=${f/SQK-RBK114-96_/}
    mv "$f" "$new_f"
done
```

Final read counts:
* O2024-029: 24126690
* O2024-030: 16065376

Due to read splitting, these counts are a bit higher than the pod5 read counts. I was originally planning on demux the pod5s as well, but due to read splitting, there isn't a one-to-one relationship, so I decided against it (also saves disk space).

Clean up:
```bash
cd /data/scratch/projects/punim1894/O2024-029
rm -r pod5
rm reads.bam

cd /data/scratch/projects/punim1894/O2024-030
rm -r pod5
rm reads.bam
```




# Read QC

https://github.com/nanoporetech/dorado/blob/release-v0.7/documentation/SAM.md

I'm just doing some basic QC here: filter on filtering on the mean basecall qscore value (in `qs:f:` tag). For my reference assemblies made with Trycycler, I may do more aggressive QC as necessary.

Get histograms of qscores:
```bash
cd /data/scratch/projects/punim1894/O2024-029
{ cat barcode02.fastq barcode03.fastq barcode04.fastq barcode05.fastq barcode06.fastq barcode07.fastq barcode09.fastq barcode13.fastq barcode14.fastq barcode15.fastq barcode16.fastq barcode25.fastq barcode26.fastq barcode28.fastq barcode29.fastq barcode30.fastq barcode31.fastq barcode32.fastq barcode41.fastq barcode50.fastq barcode51.fastq barcode52.fastq barcode56.fastq barcode57.fastq barcode63.fastq barcode68.fastq barcode69.fastq barcode75.fastq barcode77.fastq barcode78.fastq barcode79.fastq barcode80.fastq barcode82.fastq barcode83.fastq barcode84.fastq | grep -oP "qs:f:[\d\.]+" | grep -oP "[\d\.]+"; echo "0.0"; } | hist -w 0.25 -x -p '0'
{ cat unclassified.fastq | grep -oP "qs:f:[\d\.]+" | grep -oP "[\d\.]+"; echo "0.0"; } | hist -w 0.25 -x -p '0'

cd /data/scratch/projects/punim1894/O2024-030
{ cat barcode02.fastq barcode04.fastq barcode05.fastq barcode06.fastq barcode07.fastq barcode16.fastq barcode17.fastq barcode18.fastq barcode20.fastq barcode21.fastq barcode22.fastq barcode24.fastq barcode25.fastq barcode26.fastq barcode28.fastq barcode30.fastq barcode34.fastq barcode39.fastq barcode40.fastq barcode87.fastq barcode88.fastq barcode90.fastq barcode91.fastq barcode92.fastq | grep -oP "qs:f:[\d\.]+" | grep -oP "[\d\.]+"; echo "0.0"; } | hist -w 0.25 -x -p '0'
{ cat unclassified.fastq | grep -oP "qs:f:[\d\.]+" | grep -oP "[\d\.]+"; echo "0.0"; } | hist -w 0.25 -x -p '0'
```

O2024-029 barcodes with >1 Gbp of reads:
```
 384858|                                                                                             000000000                                                                                            
 364602|                                                                                           000000000000                                                                                           
 344347|                                                                                         0000000000000000                                                                                         
 324091|                                                                                        000000000000000000                                                                                        
 303835|                                                                                      000000000000000000000                                                                                       
 283580|                                                                                    000000000000000000000000                                                                                      
 263324|                                                                                   00000000000000000000000000                                                                                     
 243068|                                                                                 00000000000000000000000000000                                                                                    
 222813|                                                                               00000000000000000000000000000000                                                                                   
 202557|                                                                             00000000000000000000000000000000000                                                                                  
 182301|                                     000                                   00000000000000000000000000000000000000                                                                                 
 162046|                                    00000                                00000000000000000000000000000000000000000                                                                                
 141790|                                   0000000                            000000000000000000000000000000000000000000000                                                                               
 121534|                                  00000000                         0000000000000000000000000000000000000000000000000                                                                              
 101279|                                 00000000000                  00000000000000000000000000000000000000000000000000000000                                                                            
  81023|                                0000000000000           0000000000000000000000000000000000000000000000000000000000000000                                                                          
  60767|                               0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                        
  40512|                             0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                    
  20256|                           000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                              
      1| 0       00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 
          .   .   .   .   .   .   .   .   .   . 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 
          5   5   5   5   5   5   5   5   5   5   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   
                                                  5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   
```

O2024-029 unclassified reads:
```
 118220|                                   000                                                                                                                                                             
 111998|                                   0000                                                                                                                                                            
 105776|                                  00000                                                                                                                                                            
  99554|                                 0000000                                                                                                                                                           
  93332|                                 0000000                                                                                                                                                           
  87110|                                00000000                                                                                                                                                           
  80888|                                000000000                                                                                                                                                          
  74666|            00                 0000000000                                                                                                                                                          
  68444|            00                00000000000                                                                                                                                                          
  62222|            00               0000000000000                                                                                                                                                         
  55999|            00              00000000000000                                                                                                                                                         
  49777|            000            0000000000000000                                                                                                                                                        
  43555|           0000          000000000000000000                                                                                                                                                        
  37333|           0000        000000000000000000000                                                                                                                                                       
  31111|          000000      00000000000000000000000                                                                                                                                                      
  24889|          0000000    00000000000000000000000000                                     00000000000                                                                                                    
  18667|          000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                            
  12445|          00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                       
   6223|          000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                
      1| 0       00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 000000
        --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 
          .   .   .   .   .   .   .   .   .   . 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 
          5   5   5   5   5   5   5   5   5   5   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   
                                                  5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   
```

O2024-030 barcodes with >1 Gbp of reads:
```
 179637|                                                                                             00000000                                                                                       
 170182|                                                                                           00000000000                                                                                      
 160728|                                                                                          00000000000000                                                                                    
 151273|                                                                                        00000000000000000                                                                                   
 141819|                                                                                       0000000000000000000                                                                                  
 132364|                                                                                     0000000000000000000000                                                                                 
 122910|                                                                                    000000000000000000000000                                                                                
 113455|                                                                                  00000000000000000000000000                                                                                
 104001|                                                                                00000000000000000000000000000                                                                               
  94546|                                                                              00000000000000000000000000000000                                                                              
  85091|                                                                            00000000000000000000000000000000000                                                                             
  75637|                                                                         000000000000000000000000000000000000000                                                                            
  66182|                                                                      00000000000000000000000000000000000000000000                                                                          
  56728|                                                                    00000000000000000000000000000000000000000000000                                                                         
  47273|                                  0000                           000000000000000000000000000000000000000000000000000                                                                        
  37819|                                 00000000                    000000000000000000000000000000000000000000000000000000000                                                                      
  28364|                                00000000000             00000000000000000000000000000000000000000000000000000000000000000                                                                   
  18910|                              0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                
   9455|                            000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                          
      1| 0       000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 0
        -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 
          .   .   .   .   .   .   .   .   .   . 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 
          5   5   5   5   5   5   5   5   5   5   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 
                                                  5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5 
```

O2024-030 unclassified reads:
```
 83431|                                000                                                                                                                                                                    
 79040|                                0000                                                                                                                                                                   
 74649|                               00000                                                                                                                                                                   
 70258|                               000000                                                                                                                                                                  
 65867|                               000000                                                                                                                                                                  
 61476|                              00000000                                                                                                                                                                 
 57085|                              00000000                                                                                                                                                                 
 52694|                             0000000000                                                                                                                                                                
 48303|                             0000000000                                                                                                                                                                
 43912|                             00000000000                                                                                                                                                               
 39520|                            0000000000000                                                                                                                                                              
 35129|                            00000000000000                                                                                                                                                             
 30738|                           000000000000000                                                                                                                                                             
 26347|                          000000000000000000                                                                                                                                                           
 21956|            00           00000000000000000000                                                                                                                                                          
 17565|           000         00000000000000000000000000              00000000000000000000000000000000000000                                                                                                  
 13174|           0000       000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                             
  8783|          000000    0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                        
  4392|          0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                                
     1| 0       000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 0000 0
       ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 
         .   .   .   .   .   .   .   .   .   . 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 
         5   5   5   5   5   5   5   5   5   5   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   
                                                 5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5  
```

The qscore distribution has multiple modes:
* ~3 (mostly in unclassified reads)
* ~9 (in both classified and unclassified)
* ~24 (mostly in classified reads)

A qscore of 12.5 seems to be a good breakpoint which separates the biggest mode of the distribution (~Q24) from the rest, so I'll use that as my cutoff:
```bash
cd /data/scratch/projects/punim1894/O2024-029/reads
for b in barcode*.fastq; do
    pass=${b/.fastq/_pass.fastq}
    cat "$b" | paste - - - - | awk -F'\t' '{split($1, a, " "); for (i in a) {if (a[i] ~ /^qs:f:/) {split(a[i], b, ":"); if (b[3] >= 12.5) print $0}}}' | tr '\t' '\n' > "$pass"
done

cd /data/scratch/projects/punim1894/O2024-030/reads
for b in barcode*.fastq; do
    pass=${b/.fastq/_pass.fastq}
    cat "$b" | paste - - - - | awk -F'\t' '{split($1, a, " "); for (i in a) {if (a[i] ~ /^qs:f:/) {split(a[i], b, ":"); if (b[3] >= 12.5) print $0}}}' | tr '\t' '\n' > "$pass"
done
```



# Compress reads and clean up

To save a bit of space, I'll now gzip the reads as tightly as I can:
```bash
cd /data/scratch/projects/punim1894/O2024-029/reads
for f in *_pass.fastq; do
    sbatch --job-name=gzip --time=4:00:00 --ntasks=1 --mem=16000 --cpus-per-task=8 --wrap "pigz -11 -p8 $f"
done

cd /data/scratch/projects/punim1894/O2024-030/reads
for f in *_pass.fastq; do
    sbatch --job-name=gzip --time=2:00:00 --ntasks=1 --mem=16000 --cpus-per-task=8 --wrap "pigz -11 -p8 $f"
done
```


And I don't need the pre-QC and unclassified FASTQs (can just keep the bams):
```bash
cd /data/scratch/projects/punim1894/O2024-029/reads
rm barcode??.fastq
rm unclassified.fastq

cd /data/scratch/projects/punim1894/O2024-030/reads
rm barcode??.fastq
rm unclassified.fastq
```




# Transfer reads to Roosta

```bash
run="O2024-029"; barcode="barcode01"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096435/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode02"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097495/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode03"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097501/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode04"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097571/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode05"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097349/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode06"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097035/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode07"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096270/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode08"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097802/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode09"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097430/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode10"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00097667/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode11"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096634/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode12"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096412/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode13"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00094869/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode14"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00098114/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode15"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00095756/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode16"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096674/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode17"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00036378/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode18"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00037914/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode19"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00039769/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode20"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00023830/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode21"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00041473/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode22"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00063628/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode23"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00045646/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode25"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00021551/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode26"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027678/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode27"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00038273/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode28"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00064438/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode29"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00052630/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode30"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00028810/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode31"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00049167/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode32"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00055259/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode41"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00031899/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode42"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00018242/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode43"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00029845/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode44"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00039410/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode45"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00034285/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode46"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00037930/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode47"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00056892/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode49"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00018307/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode50"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00046200/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode51"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00017411/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode52"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010326/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode53"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00085855/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode54"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00020536/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode55"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00038884/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode56"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00030319/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode57"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00018770/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode58"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00031696/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode59"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025722/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode60"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00044976/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode61"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00049264/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode62"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00019206/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode63"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00026122/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode64"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025040/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode65"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00022702/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode66"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00044315/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode67"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00042429/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode68"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00048656/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode69"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00013571/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode70"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00035802/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode71"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00018983/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode73"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00016814/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode74"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010911/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode75"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009397/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode76"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009001/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode77"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00007207/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode78"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010001/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode79"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00007389/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode80"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009994/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode81"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009590/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode82"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009362/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode83"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009203/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-029"; barcode="barcode84"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010088/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode01"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00036386/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode02"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00026329/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode03"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00035370/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode04"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00020566/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode05"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00008243/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode06"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010052/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode07"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00028554/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode08"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00021334/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode09"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00036863/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode10"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00037801/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode11"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00008994/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode13"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00045986/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode14"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00031534/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode15"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00038111/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode16"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00018490/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode17"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010020/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode18"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025489/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode19"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025487/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode20"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00020608/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode21"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00017930/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode22"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00005735/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode23"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00029307/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode24"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00012898/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode25"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00011859/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode26"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00040532/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode27"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00037061/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode28"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00043864/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode29"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027803/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode30"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025222/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode31"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00044460/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode32"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00024150/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode33"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025114/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode34"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00029268/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode35"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00042176/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode37"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00030826/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode38"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00041863/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode39"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00091465/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode40"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00019873/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode41"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00070142/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode42"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00074515/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode43"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00004555/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode44"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00005176/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode45"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00005186/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode46"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025539/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode47"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025553/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode48"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025687/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode49"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027955/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode50"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027958/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode51"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027981/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode52"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027982/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode53"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00027954/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode54"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025527/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode55"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00008979/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode56"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025529/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode57"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00009606/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode58"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00093330/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode59"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00075502/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode61"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00040426/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode62"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00094505/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode63"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00096865/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode64"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00056693/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode65"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00058635/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode66"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00058835/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode67"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00059427/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode68"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00060270/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode69"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00065997/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode70"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00066124/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode71"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00021558/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode72"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025373/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode73"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00025378/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode74"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00023581/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode75"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00023584/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode76"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00029329/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode77"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00053859/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode78"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00054655/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode79"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00095417/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode80"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00095416/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode81"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00095408/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode82"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00072179/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode83"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00054652/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode84"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00057084/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode85"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00032793/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode86"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00031978/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode87"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00021208/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode88"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00062512/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode89"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00040126/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode90"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00015264/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode91"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00036400/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
run="O2024-030"; barcode="barcode92"; cd ~/2024-06_new_assembler_benchmark/AUSMDU00010405/reads; scp spartan:/data/scratch/projects/punim1894/"$run"/reads/"$barcode"_pass.fastq.gz nanopore.fastq.gz
```

And get some read stats:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in $(ls */reads/nanopore.fastq.gz | grep -oP "AUSMDU\d+"); do
    cd ~/2024-06_new_assembler_benchmark/"$id"/reads
    fast_count nanopore.fastq.gz | awk '{print $2 " " $3 " " $6}' | read count bases n50
    printf "$id\t$count\t$bases\t$n50\n"
done
```



# Initial sample selection

Excluded samples which met any of these criteria:
* no ONT reads (only two samples)
* Illumina reads strangely short (just one sample with half-length reads for some reason)
* too much Illumina contamination (>50kbp of low-depth sequence removed by Unicycler, only two samples)
* low ONT N50 (<5kbp, just one sample)
* low Illumina depth (<50x, just two samples)
* low ONT depth (<200x, the most common read for QC failure)

This left me with 47 samples from 16 different species.

```bash
cd ~/2024-06_new_assembler_benchmark
mkdir qc_fail
for id in AUSMDU00004555 AUSMDU00005176 AUSMDU00005186 AUSMDU00005735 AUSMDU00007389 AUSMDU00008979 AUSMDU00008994 AUSMDU00009001 AUSMDU00009590 AUSMDU00009606 AUSMDU00009994 AUSMDU00010001 AUSMDU00010020 AUSMDU00010052 AUSMDU00010405 AUSMDU00010911 AUSMDU00010915 AUSMDU00011859 AUSMDU00012898 AUSMDU00012916 AUSMDU00015264 AUSMDU00016814 AUSMDU00017930 AUSMDU00018242 AUSMDU00018307 AUSMDU00018983 AUSMDU00019206 AUSMDU00020536 AUSMDU00020608 AUSMDU00021334 AUSMDU00021558 AUSMDU00022702 AUSMDU00023581 AUSMDU00023584 AUSMDU00023830 AUSMDU00024150 AUSMDU00025040 AUSMDU00025114 AUSMDU00025373 AUSMDU00025378 AUSMDU00025487 AUSMDU00025489 AUSMDU00025527 AUSMDU00025529 AUSMDU00025539 AUSMDU00025553 AUSMDU00025687 AUSMDU00025722 AUSMDU00027803 AUSMDU00027954 AUSMDU00027955 AUSMDU00027958 AUSMDU00027981 AUSMDU00027982 AUSMDU00029307 AUSMDU00029329 AUSMDU00029845 AUSMDU00030826 AUSMDU00031534 AUSMDU00031696 AUSMDU00031978 AUSMDU00032793 AUSMDU00034285 AUSMDU00035370 AUSMDU00035802 AUSMDU00036378 AUSMDU00036386 AUSMDU00036863 AUSMDU00037061 AUSMDU00037801 AUSMDU00037914 AUSMDU00037930 AUSMDU00038111 AUSMDU00038273 AUSMDU00038884 AUSMDU00039410 AUSMDU00039769 AUSMDU00040126 AUSMDU00040426 AUSMDU00040532 AUSMDU00041473 AUSMDU00041863 AUSMDU00042176 AUSMDU00042429 AUSMDU00044315 AUSMDU00044460 AUSMDU00044976 AUSMDU00045646 AUSMDU00045986 AUSMDU00049264 AUSMDU00053859 AUSMDU00054652 AUSMDU00054655 AUSMDU00056693 AUSMDU00056892 AUSMDU00057084 AUSMDU00058635 AUSMDU00058835 AUSMDU00059427 AUSMDU00060270 AUSMDU00063628 AUSMDU00065997 AUSMDU00066124 AUSMDU00070142 AUSMDU00072179 AUSMDU00074515 AUSMDU00075502 AUSMDU00085855 AUSMDU00093330 AUSMDU00094505 AUSMDU00095408 AUSMDU00095416 AUSMDU00096412 AUSMDU00096435 AUSMDU00096634 AUSMDU00096865 AUSMDU00097667; do
    mv "$id" qc_fail
done
```




# Unicycler hybrid assembly

First, I copied over the graph from the Illumina-only assembly to save time:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    mkdir unicycler_hybrid
    cp unicycler_illumina/002_depth_filter.gfa unicycler_hybrid
done
```

I then looked at each of the `002_depth_filter.gfa` files in Bandage, checking for multiplicity problems and manually repairing them when found. Quite a lot needed fixes, both due to heterogeneity (e.g. in-and-out plasmids and fim switches) and not-great Illumina reads:
* AUSMDU00007207: no fixes
* AUSMDU00008243: fixed one long contig
* AUSMDU00009203: no fixes
* AUSMDU00009362: no fixes
* AUSMDU00009397: fixed a small plasmid with a dangling tip
* AUSMDU00010088: no fixes
* AUSMDU00010326: no fixes
* AUSMDU00013571: fixed one long contig
* AUSMDU00017411: no fixes
* AUSMDU00018490: fixed two plasmid contigs
* AUSMDU00018770: no fixes
* AUSMDU00019873: fixed two long contigs
* AUSMDU00020566: no fixes
* AUSMDU00021208: fixed a long contig with a dangling tip
* AUSMDU00021551: no fixes
* AUSMDU00025222: fixed one long contig
* AUSMDU00026122: no fixes
* AUSMDU00026329: no fixes
* AUSMDU00027678: no fixes
* AUSMDU00028554: fixed a few long contigs and plasmid contigs
* AUSMDU00028810: fixed two long contigs
* AUSMDU00029268: no fixes
* AUSMDU00030319: no fixes
* AUSMDU00031899: no fixes
* AUSMDU00036400: no fixes
* AUSMDU00043864: fixed one long contig
* AUSMDU00046200: fixed two long contigs
* AUSMDU00048656: no fixes
* AUSMDU00049167: no fixes
* AUSMDU00052630: fixed one long contig and two dangling tips
* AUSMDU00055259: no fixes
* AUSMDU00062512: fixed two long contigs
* AUSMDU00064438: fixed three long contigs
* AUSMDU00091465: no fixes
* AUSMDU00094869: no fixes
* AUSMDU00095417: no fixes
* AUSMDU00095756: no fixes
* AUSMDU00096270: no fixes
* AUSMDU00096674: no fixes
* AUSMDU00097035: no fixes
* AUSMDU00097349: no fixes
* AUSMDU00097430: no fixes
* AUSMDU00097495: no fixes
* AUSMDU00097501: no fixes
* AUSMDU00097571: no fixes
* AUSMDU00097802: no fixes
* AUSMDU00098114: fixed one long contig

Since I already did some qscore-based filtering, the only QC I'll do now is tossing out very short reads:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    cd reads_qc
    filtlong --min_length 1000 ../reads/nanopore.fastq.gz | pigz -p8 > nanopore.fastq.gz
done
```
This usually dropped about 10-20% of the reads, but since they were short, this was only 1-2% of the bases.

Launch Unicycler hybrid:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    if [ -f unicycler_hybrid/assembly.gfa ]; then continue; fi
    unicycler -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -l reads/nanopore.fastq.gz -o unicycler_hybrid -t 64
done
```

Check for complete replicons:
```bash
cd ~/2024-06_new_assembler_benchmark
for hybrid_log in AUSMDU*/unicycler_hybrid/unicycler.log; do
    id=$(echo "$hybrid_log" | sed 's|/unicycler_hybrid/unicycler.log||')
    complete_count=$(grep -cP "  complete" "$hybrid_log")
    if grep -q "incomplete" "$hybrid_log"; then
        echo "$id\t$complete_count\tno"
    else
        echo "$id\t$complete_count\tyes"
    fi
done
```




# Assembly for Trycycler

Tools:
* Canu v2.2
* Flye v2.9.4
* Miniasm v0.3
* Minipolish v0.1.3
* Racon v1.5.0
* any2fasta v0.4.2
* NECAT v20200803
* NextDenovo v2.5.2
* NextPolish v1.4.1
* Raven v1.8.3

Prepare assemblies:
```bash
declare -A genome_sizes
genome_sizes["AUSMDU00007207"]=4893123
genome_sizes["AUSMDU00008243"]=6892756
genome_sizes["AUSMDU00009203"]=5715579
genome_sizes["AUSMDU00009362"]=7046063
genome_sizes["AUSMDU00009397"]=5796608
genome_sizes["AUSMDU00010088"]=5305175
genome_sizes["AUSMDU00010326"]=5408307
genome_sizes["AUSMDU00013571"]=6097181
genome_sizes["AUSMDU00017411"]=5994004
genome_sizes["AUSMDU00018490"]=5323712
genome_sizes["AUSMDU00018770"]=6004138
genome_sizes["AUSMDU00019873"]=4391942
genome_sizes["AUSMDU00020566"]=5151652
genome_sizes["AUSMDU00021208"]=5258167
genome_sizes["AUSMDU00021551"]=5589538
genome_sizes["AUSMDU00025222"]=4902277
genome_sizes["AUSMDU00026122"]=4421204
genome_sizes["AUSMDU00026329"]=5113646
genome_sizes["AUSMDU00027678"]=4656221
genome_sizes["AUSMDU00028554"]=5148630
genome_sizes["AUSMDU00028810"]=5446907
genome_sizes["AUSMDU00029268"]=4633117
genome_sizes["AUSMDU00030319"]=6765258
genome_sizes["AUSMDU00031899"]=4961539
genome_sizes["AUSMDU00036400"]=5189554
genome_sizes["AUSMDU00043864"]=4898193
genome_sizes["AUSMDU00046200"]=6494488
genome_sizes["AUSMDU00048656"]=5834453
genome_sizes["AUSMDU00049167"]=5025502
genome_sizes["AUSMDU00052630"]=5071672
genome_sizes["AUSMDU00055259"]=5625214
genome_sizes["AUSMDU00062512"]=5334120
genome_sizes["AUSMDU00064438"]=5504511
genome_sizes["AUSMDU00091465"]=4510855
genome_sizes["AUSMDU00094869"]=3010950
genome_sizes["AUSMDU00095417"]=1664795
genome_sizes["AUSMDU00095756"]=3028725
genome_sizes["AUSMDU00096270"]=3111437
genome_sizes["AUSMDU00096674"]=2973440
genome_sizes["AUSMDU00097035"]=3082823
genome_sizes["AUSMDU00097349"]=3007252
genome_sizes["AUSMDU00097430"]=3004380
genome_sizes["AUSMDU00097495"]=3061258
genome_sizes["AUSMDU00097501"]=3029244
genome_sizes["AUSMDU00097571"]=2956679
genome_sizes["AUSMDU00097802"]=3011410
genome_sizes["AUSMDU00098114"]=3047525

cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    genome_size=$genome_sizes["$id"]
    
    cd ~/2024-06_new_assembler_benchmark/"$id"

    trycycler subsample --reads reads_qc/nanopore.fastq.gz --out_dir read_subsets --count 24 --genome_size "$genome_size"
    mkdir assemblies

    i=01
    canu -p canu -d canu_temp_"$i" -fast genomeSize="$genome_size" useGrid=false minThreads=32 maxThreads=32 maxMemory=120 -nanopore read_subsets/sample_"$i".fastq
    canu_trim.py canu_temp_"$i"/canu.contigs.fasta > assemblies/assembly_"$i".fasta
    rm -rf canu_temp_"$i"

    i=07
    canu -p canu -d canu_temp_"$i" -fast genomeSize="$genome_size" useGrid=false minThreads=32 maxThreads=32 maxMemory=120 -nanopore read_subsets/sample_"$i".fastq
    canu_trim.py canu_temp_"$i"/canu.contigs.fasta > assemblies/assembly_"$i".fasta
    rm -rf canu_temp_"$i"

    i=13
    canu -p canu -d canu_temp_"$i" -fast genomeSize="$genome_size" useGrid=false minThreads=32 maxThreads=32 maxMemory=120 -nanopore read_subsets/sample_"$i".fastq
    canu_trim.py canu_temp_"$i"/canu.contigs.fasta > assemblies/assembly_"$i".fasta
    rm -rf canu_temp_"$i"

    i=19
    canu -p canu -d canu_temp_"$i" -fast genomeSize="$genome_size" useGrid=false minThreads=32 maxThreads=32 maxMemory=120 -nanopore read_subsets/sample_"$i".fastq
    canu_trim.py canu_temp_"$i"/canu.contigs.fasta > assemblies/assembly_"$i".fasta
    rm -rf canu_temp_"$i"

    for i in 02 08 14 20; do
        flye --nano-hq read_subsets/sample_"$i".fastq --threads 32 --out-dir flye_temp
        cp flye_temp/assembly.fasta assemblies/assembly_"$i".fasta
        cp flye_temp/assembly_graph.gfa assemblies/assembly_"$i".gfa
        rm -r flye_temp
    done

    for i in 03 09 15 21; do
        miniasm_and_minipolish.sh read_subsets/sample_"$i".fastq 32 > assemblies/assembly_"$i".gfa
        any2fasta assemblies/assembly_"$i".gfa > assemblies/assembly_"$i".fasta
    done

    for i in 04 10 16 22; do
        ~/programs/NECAT/Linux-amd64/bin/necat.pl config config.txt
        realpath read_subsets/sample_"$i".fastq > read_list.txt
        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
        sed -i "s/GENOME_SIZE=/GENOME_SIZE="$genome_size"/" config.txt
        sed -i "s/THREADS=4/THREADS=12/" config.txt
        ~/programs/NECAT/Linux-amd64/bin/necat.pl bridge config.txt
        cp necat/6-bridge_contigs/polished_contigs.fasta assemblies/assembly_"$i".fasta
        rm -r necat config.txt read_list.txt
    done

    for i in 05 11 17 23; do
        echo read_subsets/sample_"$i".fastq > input.fofn
        cp ~/programs/NextDenovo/doc/run.cfg nextdenovo_run.cfg
        sed -i "s/genome_size = 1g/genome_size = "$genome_size"/" nextdenovo_run.cfg
        sed -i "s/parallel_jobs = 20/parallel_jobs = 1/" nextdenovo_run.cfg
        sed -i "s/read_type = clr/read_type = ont/" nextdenovo_run.cfg
        sed -i "s/pa_correction = 3/pa_correction = 1/" nextdenovo_run.cfg
        sed -i "s/correction_options = -p 15/correction_options = -p 32/" nextdenovo_run.cfg
        sed -i "s/-t 8/-t 32/" nextdenovo_run.cfg
        nextDenovo nextdenovo_run.cfg
        cp 01_rundir/03.ctg_graph/nd.asm.fasta nextdenovo_temp.fasta
        rm -r 01_rundir nextdenovo_run.cfg input.fofn
        echo read_subsets/sample_"$i".fastq > lgs.fofn
        cat ~/programs/NextPolish/doc/run.cfg | grep -v "sgs" | grep -v "hifi" > nextpolish_run.cfg
        sed -i "s/parallel_jobs = 6/parallel_jobs = 1/" nextpolish_run.cfg
        sed -i "s/multithread_jobs = 5/multithread_jobs = 32/" nextpolish_run.cfg
        sed -i "s|genome = ./raw.genome.fasta|genome = nextdenovo_temp.fasta|" nextpolish_run.cfg
        sed -i "s|-x map-ont|-x map-ont -t 32|" nextpolish_run.cfg
        nextPolish nextpolish_run.cfg
        cp 01_rundir/genome.nextpolish.fasta assemblies/assembly_"$i".fasta
        rm -r 01_rundir pid*.log.info nextpolish_run.cfg lgs.fofn nextdenovo_temp.fasta
    done

    for i in 06 12 18 24; do
        raven --threads 32 --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_"$i".gfa read_subsets/sample_"$i".fastq > assemblies/assembly_"$i".fasta
    done
done
```



# Trycycler

Trycycler cluster:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/nanopore.fastq.gz --out_dir trycycler --threads 64
done
```

And then Trycycler reconcile:
```bash
trycycler reconcile --reads reads_qc/nanopore.fastq.gz --cluster_dir trycycler/cluster_xxx  # run for each good cluster
```

And the remaining Trycycler steps:
```bash
for c in trycycler/cluster_*; do trycycler msa --threads 64 --cluster_dir "$c"; done
trycycler partition --reads reads_qc/nanopore.fastq.gz --cluster_dirs trycycler/cluster_* --threads 96
for c in trycycler/cluster_*; do trycycler consensus --cluster_dir "$c"; done
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta
rm trycycler/cluster_*/4_reads.fastq
```













# Illumina polishing


Tools:
* Polypolish v0.6.0
* Pypolca v0.3.1

```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    echo "\n\n\n\n\n\n\n$id"
    cd ~/2024-06_new_assembler_benchmark/"$id"
    bwa index trycycler.fasta
    bwa mem -t 24 -a trycycler.fasta reads_qc/illumina_1.fastq.gz > alignments_1.sam
    bwa mem -t 24 -a trycycler.fasta reads_qc/illumina_2.fastq.gz > alignments_2.sam
    polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish polish trycycler.fasta filtered_1.sam filtered_2.sam > trycycler_polypolish.fasta
    rm *.amb *.ann *.bwt *.pac *.sa *.sam
    pypolca run --careful -a trycycler_polypolish.fasta -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -t 24 -o pypolca
    cp pypolca/pypolca_corrected.fasta trycycler_polypolish_pypolca.fasta
    rm -r pypolca
done
```

Count differences:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    compare_assemblies.py --aligner edlib trycycler.fasta trycycler_polypolish.fasta > polypolish.changes
    compare_assemblies.py --aligner edlib trycycler_polypolish.fasta trycycler_polypolish_pypolca.fasta > pypolca.changes
    compare_assemblies.py --aligner edlib trycycler.fasta trycycler_polypolish_pypolca.fasta > all_polishing.changes
done

for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    polypolish=$(grep -o "*" polypolish.changes | wc -l)
    pypolca=$(grep -o "*" pypolca.changes | wc -l)
    all_polishing=$(grep -o "*" all_polishing.changes | wc -l)
    size=$(fast_count trycycler_polypolish_pypolca.fasta | cut -f3)
    printf "$polypolish\t$pypolca\t$all_polishing\t$size\n"
done
```

This all seemed to polish up nicely, so I'm happy to keep the results as a ground-truth reference:
```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    seqtk seq trycycler_polypolish_pypolca.fasta > reference.fasta
    ~/2024-06_new_assembler_benchmark/scripts/rename_contigs.py reference.fasta
done
```



# Mash tree

This is to check for redundancy in the genomes.

```bash
cd ~/2024-06_new_assembler_benchmark

threads=12
sketch_size=10000

mash sketch -p $threads -o mash_sketches -s $sketch_size AUS*/reference.fasta
mash dist -p $threads mash_sketches.msh mash_sketches.msh -t > distances.tab
tail -n +2 distances.tab > distances.tab.temp  # remove first line
wc -l distances.tab.temp | awk '"[0-9]+ errors" {sum += $1}; END {print sum}' > distances.ndist
cat distances.ndist distances.tab.temp > mash.phylip
rm distances.ndist mash_sketches.msh distances.tab distances.tab.temp
```

```r
library(ape)
library(phangorn)
distances <- readDist('mash.phylip')
tree <- fastme.bal(distances)
tree$edge.length <- pmax(tree$edge.length, 0.0)  # set any negative branch lengths to zero
tree <- midpoint(tree)
write.tree(tree, 'mash.newick')
```

The only two genomes which are extremely close are AUSMDU00019873 and AUSMDU00091465 (Mash distance of 0.000264). So I dropped AUSMDU00019873 because it had lower ONT depth and more Illumina polishing changes.




# Check for heterogeneity

Tools:
* Sniffles v2.4

```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    minimap2 -a -x map-ont -t 48 reference.fasta reads_qc/nanopore.fastq.gz | samtools sort > nanopore.bam
    samtools index nanopore.bam
    sniffles -i nanopore.bam -v sniffles.vcf
    rm nanopore.bam*
done
```

Check to see which samples had a result:
```bash
cd ~/2024-06_new_assembler_benchmark
for v in AUSMDU*/sniffles.vcf; do
    echo "\n\n$v"
    cat "$v" | grep "AF=" | grep -vP "\t1\tSniffles2"
done
```

Eight samples had a Sniffles variant, so I excluded those, leaving me with 18 samples.






# Read cleanup

In this step, I ensure that each read aligns to the reference. This mainly serves to clean up reads which are total junk or cross-barcode contamination. For example, some isolates in these runs had very-high-depth small plasmids, which could potentially be present in other barcodes at assemblable depths.

To pass the filter, at least half of the read's length needs to align to some part of the reference with any identity.

```bash
cd ~/2024-06_new_assembler_benchmark
for id in AUSMDU*; do
    cd ~/2024-06_new_assembler_benchmark/"$id"
    cat reference.fasta | awk '/^>/ {print} !/^>/ {print $0$0}' > doubled_reference.fasta
    minimap2 -c -x map-ont -t 48 doubled_reference.fasta reads/nanopore.fastq.gz > doubled_reference.paf
    ~/2024-06_new_assembler_benchmark/scripts/read_filter.py reads/nanopore.fastq.gz doubled_reference.paf | pigz -p8 > all_reads.fastq.gz
    rm doubled_reference.*
done
```

This resulted in up to ~2% of the reads being dropped, but they must have been mostly short reads, as only the number of bases only went down by ~0.1%.





# Done!

In the end, I have these 18 genomes:
* AUSMDU00009397 Shigella flexneri
* AUSMDU00017411 Klebsiella pneumoniae
* AUSMDU00018770 Klebsiella pneumoniae
* AUSMDU00021551 Enterobacter hormaechei
* AUSMDU00026122 Providencia rettgeri
* AUSMDU00027678 Enterobacter hormaechei
* AUSMDU00049167 Enterobacter kobei
* AUSMDU00091465 Salmonella enterica
* AUSMDU00094869 Listeria monocytogenes
* AUSMDU00095417 Campylobacter jejuni
* AUSMDU00096674 Listeria monocytogenes
* AUSMDU00097349 Listeria welshimeri
* AUSMDU00097430 Listeria innocua
* AUSMDU00097495 Listeria monocytogenes
* AUSMDU00097501 Listeria monocytogenes
* AUSMDU00097571 Listeria seeligeri
* AUSMDU00097802 Listeria monocytogenes
* AUSMDU00098114 Listeria monocytogenes

Each genomes has:
* A good reference with no ambiguity
  * Trycycler assembly with no complications
  * Illumina polishing which made no more than 4 changes
  * In good agreement with a Unicycler hybrid assembly
  * No structural variants found by Sniffles
* A good read set
  * At least 200x depth
  * No junk or contaminant reads
