cd ../build

trajectories=$1
number_threads=${2:-1}
trajectories_per_thread=`expr $trajectories / $number_threads`

echo "Starting $number_threads threads..."

for thread in {1..$number_threads}
do
    seed=$((41 + 100 * $thread))
    ./otg-test $trajectories_per_thread $seed --abort-after=5 -nv &
done

wait
