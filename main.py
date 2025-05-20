from fraction import Fraction

def read_from_file(filename: str = "input.txt"):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    z_line = lines[0].split()
    z = [Fraction(int(coeff)) for coeff in z_line]

    matrix = []
    for line in lines[1:]:
        row = [Fraction(int(num)) for num in line.split()]
        matrix.append(row)
    return {'z': z, 'matrix': matrix}

def get_basic_vars(matrix):
    basic_vars = []
    num_rows = len(matrix)
    num_cols = len(matrix[0]) - 1
    for col in range(num_cols):
        ones = [i for i in range(num_rows) if matrix[i][col] == Fraction(1)]
        zeros = [i for i in range(num_rows) if matrix[i][col] == Fraction(0)]
        if len(ones) == 1 and len(zeros) == num_rows - 1:
            basic_vars.append(col)
    return basic_vars


def print_simplex_table(simplex_table, basic_vars, co_row=None):
    num_rows = len(simplex_table) - 1
    num_cols = len(simplex_table[0])
    field_width = 10
    total_width = 7 + (field_width + 1) * num_cols

    print("=" * total_width)
    print("{:^6}|".format("б.п."), end='')
    print("{:^{size}}|".format("1", size=field_width), end='')
    for i in range(num_cols - 1):
        print("{:^{size}}|".format("x" + str(i + 1), size=field_width), end='')
    print()
    print("-" * total_width)

    for i in range(num_rows):
        print("{:^6}|".format("x" + str(basic_vars[i] + 1)), end='')
        for j in range(num_cols):
            value = str(simplex_table[i][j])
            print("{:^{size}}|".format(value, size=field_width), end='')
        print()

    print("-" * total_width)
    print("{:^6}|".format("Z"), end='')
    for j in range(num_cols):
        value = str(simplex_table[-1][j])
        print("{:^{size}}|".format(value, size=field_width), end='')
    print()

    # Вычисление строки CO, если она не передана
    if co_row is None:
        b = [simplex_table[i][0] for i in range(num_rows)]
        negative_b = [i for i, val in enumerate(b) if val < Fraction(0)]
        if negative_b:
            pivot_row = min(negative_b, key=lambda i: b[i])
            co_row = []
            for j in range(1, num_cols):
                if simplex_table[pivot_row][j] < Fraction(0):
                    ratio = abs(simplex_table[-1][j] / simplex_table[pivot_row][j])
                    co_row.append(ratio)
                else:
                    co_row.append("-")
        else:
            co_row = ["-"] * (num_cols - 1)

    print("-" * total_width)
    print("{:^6}|".format("CO"), end='')
    print("{:^{size}}|".format("-", size=field_width), end='')
    for j in range(1, num_cols):
        value = str(co_row[j - 1]) if co_row[j - 1] != "-" else "-"
        print("{:^{size}}|".format(value, size=field_width), end='')
    print()

    print("=" * total_width)


def dual_simplex_method(simplex_table, basic_vars):
    num_rows = len(simplex_table) - 1
    num_cols = len(simplex_table[0])

    if any(simplex_table[-1][j] < Fraction(0) for j in range(1, num_cols)):
        print("\nНачальная симплекс-таблица:\n")
        print_simplex_table(simplex_table, basic_vars)
        print("\nУсловие двойственности не выполнено: есть отрицательные коэффициенты в Z-строке.")
        print("Двойственный симплекс-метод не применим.")
        return

    print("\nНачальная симплекс-таблица:\n")
    b = [simplex_table[i][0] for i in range(num_rows)]
    negative_b = [i for i, val in enumerate(b) if val < Fraction(0)]
    if negative_b:
        pivot_row = min(negative_b, key=lambda i: b[i])
        co_row = []
        for j in range(1, num_cols):
            if simplex_table[pivot_row][j] < Fraction(0):
                ratio = simplex_table[-1][j] / (-simplex_table[pivot_row][j])
                co_row.append(ratio)
            else:
                co_row.append("-")
        print_simplex_table(simplex_table, basic_vars, co_row=co_row)
    else:
        print_simplex_table(simplex_table, basic_vars)

    iteration = 0
    while True:
        iteration += 1
        print(f"\nИтерация {iteration}:\n")

        b = [simplex_table[i][0] for i in range(num_rows)]
        negative_b = [i for i, val in enumerate(b) if val < Fraction(0)]
        if not negative_b:
            print_simplex_table(simplex_table, basic_vars)
            print("Оптимальное решение найдено.")

            non_basis_zero = [j for j in range(1, num_cols) if
                              (j - 1) not in basic_vars and simplex_table[-1][j] == Fraction(0)]
            if non_basis_zero:
                print("Обнаружено множество решений.")
                entering_var = non_basis_zero[0] # Находим второе базисное решение
                ratios = []
                for i in range(num_rows):
                    if simplex_table[i][entering_var] > Fraction(0):
                        ratio = simplex_table[i][0] / simplex_table[i][entering_var]
                        ratios.append((ratio, i))
                if not ratios:
                    print("Не удалось найти второе базисное решение.")
                    return
                ratios.sort()
                leaving_row = ratios[0][1]
                new_table = [row.copy() for row in simplex_table]
                pivot_val = new_table[leaving_row][entering_var]
                for j in range(num_cols):
                    new_table[leaving_row][j] /= pivot_val
                for i in range(num_rows + 1):
                    if i == leaving_row:
                        continue
                    factor = new_table[i][entering_var]
                    for j in range(num_cols):
                        new_table[i][j] -= factor * new_table[leaving_row][j]
                new_basic_vars = basic_vars.copy()
                new_basic_vars[leaving_row] = entering_var - 1
                print("\nСимплекс таблица для второго базисного решения:")
                print_simplex_table(new_table, new_basic_vars)
                solution1 = {basic_vars[i]: simplex_table[i][0] for i in range(num_rows)}
                solution2 = {new_basic_vars[i]: new_table[i][0] for i in range(num_rows)}
                print("Первое базисное решение:")
                for var in range(num_cols - 1):
                    if var in solution1:
                        print(f"x{var + 1} = {solution1[var]}")
                    else:
                        print(f"x{var + 1} = 0")
                print("Второе базисное решение:")
                for var in range(num_cols - 1):
                    if var in solution2:
                        print(f"x{var + 1} = {solution2[var]}")
                    else:
                        print(f"x{var + 1} = 0")
                return
            else:
                solution = {basic_vars[i]: simplex_table[i][0] for i in range(num_rows)}
                for var in range(num_cols - 1):
                    if var in solution:
                        print(f"x{var + 1} = {solution[var]}")
                    else:
                        print(f"x{var + 1} = 0")
                print(f"Значение Z_max = {simplex_table[-1][0]}")
                return

        pivot_row = min(negative_b, key=lambda i: b[i])

        co_row = []
        for j in range(1, num_cols):
            if simplex_table[pivot_row][j] < Fraction(0):
                ratio = simplex_table[-1][j] / (-simplex_table[pivot_row][j])
                co_row.append(ratio)
            else:
                co_row.append("-")

        if all(simplex_table[pivot_row][j] >= Fraction(0) for j in range(1, num_cols)):
            print_simplex_table(simplex_table, basic_vars, co_row=co_row)
            print("\nНет решений: в разрешающей строке нет отрицательных коэффициентов.")
            return

        ratios = [(co_row[j - 1], j) for j in range(1, num_cols) if co_row[j - 1] != "-"]
        if not ratios:
            print("\nНе удалось найти разрешающий столбец.")
            return
        min_ratio, pivot_col = min(ratios)

        print_simplex_table(simplex_table, basic_vars, co_row=co_row)
        print(f"Разрешающая строка: x{basic_vars[pivot_row] + 1} (b = {simplex_table[pivot_row][0]})")
        print(f"Разрешающий столбец: x{pivot_col} (минимальное отношение = {min_ratio})")
        print(f"Разрешающий элемент: {simplex_table[pivot_row][pivot_col]}")

        pivot_val = simplex_table[pivot_row][pivot_col]
        for j in range(num_cols):
            simplex_table[pivot_row][j] /= pivot_val
        for i in range(num_rows + 1):
            if i == pivot_row:
                continue
            factor = simplex_table[i][pivot_col]
            for j in range(num_cols): # правило квадрата упрощено - просто делим строку
                simplex_table[i][j] -= factor * simplex_table[pivot_row][j]
        basic_vars[pivot_row] = pivot_col - 1

if __name__ == '__main__':
    """Основная функция."""
    data = read_from_file()
    z, matrix = data["z"], data["matrix"]
    basic_vars = get_basic_vars(matrix)
    simplex_table = [[row[-1]] + row[:-1] for row in matrix] + [[z[-1]] + z[:-1]]
    dual_simplex_method(simplex_table, basic_vars)

