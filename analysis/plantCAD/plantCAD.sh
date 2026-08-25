python src/zero_shot_score.py \
    -input-vcf examples/PlantAD_Rice.vcf \
    -input-fasta examples/HuangHuaZhan.fixed.fa \
    -output examples/RICE_scored_variants_cad2.vcf \
    -model 'kuleshov-group/PlantCAD2-Large-l48-d1536' \
    -contextSize 8192 \
    -batchSize 2 \
    -device 'cuda:0'