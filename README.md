# mdkp
Для корректной работы необходимы [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/) и [Boost C++ Libraries](https://www.boost.org).
Команды нужно выполнять из корневого каталога репозитория.

Собрать исполняемый файл:
```
rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . && cd ..
```
Запустить тесты `benchmarks/70x30_auto`, имя которых соотвествуюет регулярному выражению `\d+$`, на алгоритм `coral` с `8` потоками и ограничением на размер ядра `10`:
```
python3 scripts/run_tests.py -c './build/mdkp --solver_type coral --thread_count 8 --core_size 10' -p 'benchmarks/70x30_auto' -f '\d+$'
```
