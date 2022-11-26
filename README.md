# Реализация метода вычисления корней кубического многочлена с помощью аналитической формулы из статьи "Analytical formula for the roots of the general complex cubic polynomial" Автора Ibrahim Baydoun 

Представлены две реализации метода.

1. Полный повтор статьи
2. С использованием функции fma для более точных вычислений

# Эксперименты

Прогнал алгоритм 10 раз по 1'000'000 примеров, получил следующие абсолютные погрешности:

1. max absolute error 0.598953
2. max absolute error 0.598954
3. max absolute error 0.598953
4. max absolute error 0.59895
5. max absolute error 2.46281 ??? 
6. max absolute error 0.598952
7. max absolute error 0.598951
8. max absolute error 0.598952
9. max absolute error 0.598952
10. max absolute error 0.598953

# Далее

- Из всех экспериментов выделился 5 - ый с погрешностью 2.46281. Нужно обнаружить корни, которые 'неудобны' для данного алгоритма 
- Попытаться уменьшить величину погрешности
