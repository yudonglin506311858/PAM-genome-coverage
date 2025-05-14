# PAM-genome-coverage

#软件功能：给定任意基因组序列和PAM，可以计算该PAM结合到基因组的位置，数量和比例。

#输入：
-genome 具体的基因组序列文件
-pam 具体的PAM如NGG
-t 具体使用的线程

#输出：
-o 输出的文件名如results_demo.txt
tail  results_demo.txt查看具体的比例


#示例：

cd /data/yudonglin/software/genome_coverage

python /data/yudonglin/software/genome_coverage/genome_coverage_multicore.py  -genome demo.fa -pam NGG -o ngg_results_demo.txt  -t 60
![image](https://github.com/user-attachments/assets/a67e2531-6e1e-4d0b-a7a0-18b8ead9ac4f)

python /data/yudonglin/software/genome_coverage/genome_coverage_multicore.py  -genome demo.fa -pam A -o A_results_demo.txt  -t 60
![image](https://github.com/user-attachments/assets/18263797-34ca-4f5d-a0f3-03461855d316)



