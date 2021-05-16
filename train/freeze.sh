for i in `seq 4`
do
cd $i
python -m deepmd freeze 
mv frozen_graph.pb graph${i}.pb
cd ..
done

