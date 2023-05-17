package ots_calc
//import org.kotlinmath.Complex
//import org.kotlinmath.R
//import org.kotlinmath.complex
import org.kotlinmath.*
import java.lang.Exception
import java.util.*

/**
 * Класс производящий расчет токов и напряжений
 *
 * @param tracks - массив путей
 * @param meshes - массив сеток
 * @property
 */
class Calc(
    private val tracks: Array<Track>,
    private val meshes: Array<Mesh>,
    private val mpss: Array<Mps>,
)
{

    private lateinit var I_mps: Array<Real> // массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
    private lateinit var num_track_mps: Array<IntArray> // двумерный номеров путей для МПС (в каждой строке начальная и конечная точка)
    private lateinit var index_mps: Array<IntArray> //двумерный массивы  индексов узлов по сетке МПС (в каждой строке начальная и конечная точка)
    private lateinit var a_x_find: Array<Array<Real>> // матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя.
    private lateinit var U_const: Array<Real> //  массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения
    val computing_settings: Computing_settings = Computing_settings()
    val err_and_mes: Errors_and_messages = Errors_and_messages()
    init{
        verify_data(true, true, true)
        tracks.forEach { it.m3db = it.mesh.create_3diag_matrix_band(it) }
        init_find_index_mps()
        init_eval_a_x_find()
    }
    /**
     * Фунция проводит набор проверок корректности исходных данных
     */
    fun check(){
        println("Refactoring")
    }

    /**
     * TODO("add fun definition")
     */
    private fun init_find_index_mps(){
        num_track_mps = Array(mpss.size) { IntArray(2) } // инициализируем массивы
        index_mps = Array(mpss.size) { IntArray(2) }
        a_x_find =Array(mpss.size) { Array<Real>(mpss.size){0.R} } // инициализируем общий массив влияния
        for (i in 0 until mpss.size) { //заполняем массивы
            num_track_mps[i][0] = mpss[i].startTrack
            num_track_mps[i][1] = mpss[i].endTrack
            index_mps[i][0] = tracks[mpss[i].startTrack].mesh.find_near_index_over_mesh(mpss[i].startPoint)
            index_mps[i][1] = tracks[mpss[i].endTrack].mesh.find_near_index_over_mesh(mpss[i].endPoint)
        }

    }
    /** Метод для решения СЛАУ с 3диагонал ленточной матрицей
     * методом двойной проходки
     * на входе матрица и вектор правой части
     * на выходе вектор ответов
     */
    private fun solve_3diag_band(matrixBand: Array<Array<Real>>, vectorB: Array<Real> ): Array<Real> { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части
        val N = vectorB.size
        val v = Array<Real>(N){0.R}
        val u = Array<Real>(N){0.R}
        val out = Array<Real>(N){0.R}
        //прямая проходка
        v[0] = matrixBand[2][1] / -matrixBand[1][0]
        u[0] = -vectorB[0] / -matrixBand[1][0]
        for (i in 1 until N - 1) {
            v[i] = matrixBand[2][i + 1] / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
            u[i] = (matrixBand[0][i - 1] * u[i - 1] - vectorB[i]) / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
        }
        v[N - 1] = 0.0.R
        u[N - 1] = (matrixBand[0][N - 2] * u[N - 2] - vectorB[N - 1]) / (-matrixBand[1][N - 1] - matrixBand[0][N - 2] * v[N - 2])
        //обратная проходка
        out[N - 1] = u[N - 1]
        for (i in N - 1 downTo 1) {
            out[i - 1] = out[i] * v[i - 1] + u[i - 1]
        }
        return out
    }
    /**
     * Фукция расчета коэффициентов влиияния МПС от них самих
     */
    private fun init_eval_a_x_find(){
        val N_find = this.mpss.size /* количестов поисковых точек */
        val Imps = 1000.0.R
        // уловный ток МПС для определения коэффициентов влияния
        var u1: Array<Real>
        var u2: Array<Real>
        // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
        val a1_x_find = Array(N_find) { Array<Real>(N_find){0.R} }
        // временные массивы влияния в МПС в начальной и
        val a2_x_find = Array(N_find) { Array<Real>(N_find){0.R} } // в конечной
        // Перебераем поисковые точки МПС
        for (i in 0 until N_find) {
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] = -Imps // задаём ток в начальной точке подключения в данном пути
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] = Imps // задаём ток в конечной точке подключения в данном пути

            u1 = solve_3diag_band(
                tracks[num_track_mps[i][0]].m3db,
                tracks[num_track_mps[i][0]].vectorB
            ) //снимаем напряжение по сетке в пути для начальной точки подключения
            u2 = solve_3diag_band(
                tracks[num_track_mps[i][1]].m3db,
                tracks[num_track_mps[i][1]].vectorB
            ) //снимаем напряжение по сетке в пути для конечной точки подключения
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] = 0.R // обнуляем ток в текущей начальной точке подключения
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] = 0.R // обнуляем ток в текущей конечной точке подключения
            for (j in 0 until N_find) { //снова проходим по всем точкам МПС
                if (num_track_mps[i][0] == num_track_mps[j][0]) { // если номер пути начальной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u1[index_mps[j][0]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][1] == num_track_mps[j][0]) { // если номер пути конечной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u2[index_mps[j][0]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][0] == num_track_mps[j][1]) { // если номер пути начальной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u1[index_mps[j][1]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][1] == num_track_mps[j][1]) { // если номер пути конечной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u2[index_mps[j][1]] / Imps // находим коэффициенты влияния для текущей точки
                }
                a_x_find[i][j] = a1_x_find[i][j] - a2_x_find[i][j]
            }
        }
    }
    /**
     * Заполняет вектор правой части по заданному двумерному массиву с коэффициентом тока для точек с одним узлом
     * @param mesh - сетка
     * @param vectorB - вектор правой части
     * @param arr_XI - массив PV
     * @param coeff_I - коэффициент на который умножается значение тока
     */
    private fun addValues_vector_b_1node( mesh: Mesh, vectorB: Array<Real>, arr_XI: Array<PV>, coeff_I: Real ): Array<Real>
    {
        for (i in arr_XI.indices) {
            vectorB[mesh.find_near_index_over_mesh(arr_XI[i].point)] += arr_XI[i].value * coeff_I
        }
        return vectorB
    }
    /** Добавляет значения в вектор правой части по заданному двумерному массиву с коэффициентом тока для точек с двумя узлами
     * @param mesh - сетка
     * @param vectorB - вектор правой части А
     * @param arr_XI -  массив PV
     * @param coeff_I  - коэффициент тока
     */
    private fun addValues_vector_b_2node( mesh: Mesh, vectorB: Array<Real>, arr_XI: Array<PV>, coeff_I: Real  ): Array<Real>
    {
        for (i in arr_XI.indices) {
            val index2 = mesh.find_2near_index_over_mesh( arr_XI[i].point ) // номер левого и правого узла
            val k = (arr_XI[i].point - mesh.X[0]) / mesh.dX // дробный индекс точки отнсительно номеров узлов сетки
            val I2 = arr_XI[i].value * coeff_I * (k - index2[0]) // ток в левый узел
            val I1 = arr_XI[i].value - I2 // ток в правый узел
            vectorB[index2[0]] += I1
            vectorB[index2[1]] += I2
        }
        return vectorB
    }
    /**
     *  Рассчит напряжение в узлах каждого пути от постоянных источников
     */
    private fun eval_U_const() {
        for (tr in tracks){
            tr.vectorB = addValues_vector_b_1node(tr.mesh, tr.vectorB, tr.fot, -1.0.R ) // в точках ФОТ ток в одном узле с минусом
            tr.vectorB = addValues_vector_b_2node(tr.mesh, tr.vectorB, tr.eps, 1.0.R ) // в точках ЭПС ток  в двух ближайших узлах
            tr.U = solve_3diag_band(tr.m3db, tr.vectorB) //потенциал в рельсах в узлах сетки
            println("tr.name.u, ${Arrays.toString(tr.U) } ")
        }
        U_const = Array<Real>(mpss.size){ j -> tracks[num_track_mps[j][0]].U[index_mps[j][0]] - tracks[num_track_mps[j][1]].U[index_mps[j][1]]}
        println("U_const:, ${Arrays.deepToString(U_const)} ")
    }

    /**
     * Заполняет вектор правой части всех точек нулями
     */
    private fun zeros_vector_b() {
        // проход по всем путям обнуление правой части ФОТ и ЭПС
        for ( tr in tracks ){
            // обнуление для точек в один узел сетки ФОТ
            for (fotPoint in tr.fot ) {
                tr.vectorB[tr.mesh.find_near_index_over_mesh(fotPoint.point)] = 0.0.R
            }
            //обнуление для точек в два узла сетки ЭПС
            for (e in tr.eps) {
                val inds = tr.mesh.find_2near_index_over_mesh(e.point)
                tr.vectorB[inds[0]] = 0.0.R
                tr.vectorB[inds[1]] = 0.0.R
            }
        }
        // проход по всем МПС
        for (i in mpss.indices) {
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] = 0.0.R
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] = 0.0.R
        }
    }
    /**
     * Метод для расчёта мгновенной схемы ОТС с указанием начального значения тока в МПС
     */
    fun calc_ots(init_I_poisk: Array<Real>): Boolean {
        return if (!verifi_I_poisk(init_I_poisk)) {
            false
        } else calc_I_poisk(init_I_poisk)
    }

    /**
     * Метод для расчёта мгновенной схемы ОТС без указания начального значения тока в поисковых точках
     */
    fun calc_ots(): Boolean {
        val init_I_poisk = Array<Real>(mpss.size){0.R}
        return calc_I_poisk(init_I_poisk)
    }
    /**
     * По сути рассчитывает токи в поисковых точках: МПС. После этого все граничные условия во всех точках с втекающим током определены
     * Расчёт каждого поискового тока сводится к вычислению невязки напряжения на данном элементе и по величине невязки корректируется ток элемента
     * По сути рассчитывается алгебраическая система уравнений методом Ньютона в цикле итераций
     */
    private fun calc_I_poisk(init_I_poisk: Array<Real>): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        err_and_mes.reset_solver_error() //обнулим ошибки решателя
        err_and_mes.calc_completed = true
        verify_data( false, false, true ) // проверка исходных данных только координ ЭПС, т.к. остальное проверено при инициализации

        if (err_and_mes.data_error) { //проверка если данные корректны
            err_and_mes.solver_error = true
            err_and_mes.messeg_solver_error = "Расчёт невозможен. Ошибка исходных данных"
            return false
        }
        val N_mps = mpss.size //количество МПС
        val U_find = Array<Real>(N_mps){0.R} // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
        val resid_U_mps = Array<Real>(N_mps){0.R}
        var mean_resid: Double             // средняя невязка
        val limit_mean_resid = computing_settings.convergence_U     // задаём предельную невязку по достижении которой сходимость из класса computing_settings
        var damping_factor = computing_settings.initial_damping_factor        //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
        var mean_resid_pred: Double        /* значение невязки по напряжению на предыдущем шаге итераций */
        var iter = 0
        val iter_max = computing_settings.max_number_iterat // счётчик итераций и максимальное число итераций
        var counter_not_exceeded: Boolean
        var convergence_not_achieved: Boolean // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций
        eval_U_const() // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counter_not_exceeded = true
        convergence_not_achieved = true
        mean_resid = 100.0 //начальное значение средняя невязка до первой итерации
        while (counter_not_exceeded && convergence_not_achieved) {
            mean_resid_pred = mean_resid //предыдущая невязка обновление
            mean_resid = 0.0 //текущая невязка скидывается
            for (i in 0 until N_mps) {
                U_find[i] = U_const[i] //начинаем с постоянного напряжения (от заданных источников тока ФОТ ЭПС)
                for (j in 0 until N_mps) {
                    U_find[i] += init_I_poisk[j] * a_x_find[i][j] //добавляем напряжение от МПС
                }
                resid_U_mps[i] = U_find[i] - init_I_poisk[i] * mpss[i].resValue//невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                init_I_poisk[i] += damping_factor * resid_U_mps[i] / (-0.5 * a_x_find[i][i] + mpss[i].resValue) //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования
                mean_resid += resid_U_mps[i].mod //обновляем невязку
            }
            mean_resid = mean_resid / N_mps //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (mean_resid > mean_resid_pred  ) { //FIX ME - тут как то по умному нужно считать невязку в комплексной плоскости
                    damping_factor = damping_factor * 0.7
                }
            }
            iter += 1 //обновляем счётчик итераций
            counter_not_exceeded = iter < iter_max  // обновляем булевые переменные выхода из цикла итераций
            convergence_not_achieved = mean_resid > limit_mean_resid
            //println("iter=$iter mean_resid=$mean_resid damping_factor=$damping_factor")
        }
        computing_settings.current_state_solver = doubleArrayOf( (iter - 1).toDouble(), mean_resid, damping_factor ) // записываем текущее состояние решателя
        println("I_mps:, ${Arrays.toString(init_I_poisk)} ")
        I_mps = init_I_poisk // заносим токи в МПС в массивы родительского класса
        eval_node_from_all_I()
        zeros_vector_b()
        if (convergence_not_achieved) {
            err_and_mes.solver_error = true
            err_and_mes.messeg_solver_error = "Превышено максимальное число итераций равное " + computing_settings.max_number_iterat + " заданная сходимость (сред невязка по напряжен) равная " + computing_settings.convergence_U + " не достигнута."
        }
        return !convergence_not_achieved // возврат обратное к несходимости: истина если расчёт сошёлся
    }
    /**
     * При найденных токах МПС рассчитывает напряжения и ток в рельсах решением М3ДЛ
     */
    private fun eval_node_from_all_I() {
        var g_rh: Real
        var g_lf: Real //условная проводимость слева и справа от узла сетки на схеме дискретизации рельсов
        // FIX ME it never usedvar num_mesh: Int  // номер сетки для начального и конечного пути МПС
        var N: Int // количество узлов сетки

        // добавляем в вектор правой часть токи МПС (токи ФОТ и ЭПС добавлены процедурой eval_U_const() )
        for (i in mpss.indices) {
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] -= I_mps[i]
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] += I_mps[i]
        }

        //потенциал в рельсах в узлах сетки для каждого пути
        for (tr in tracks) {
            tr.U = solve_3diag_band(tr.m3db, tr.vectorB)
        }

        //расчёт токов в рельсах и тока в земле в узлах сетки для каждого пути
        for (tr in tracks) {
            N = tr.mesh.mesh_N // число узлов сетки для текущего пути
            //ток в рельсах и земле для первого узла
            g_lf = 1 / tr.Rv0 //проводимость слева от первого узла через волновое сопротивление в начале
            g_rh = -1 * tr.m3db[0][0] //проводимость справа от первого узла через нижнюю диагональ первый элемент
            tr.I[0] = 0.5 * ((0 - tr.U[0]) * g_lf + (tr.U[0] - tr.U[1]) * g_rh) //ток первого узла как полусумма токов слева и справа от него
            tr.Ignd[0] = tr.U[0] * (tr.m3db[1][0] + tr.m3db[0][0]) //ток в земле для первого узла
            //ток в рельсах и земле для остальных узлов со второго до предпоследнего
            for (i in 1 until N - 1) {
                g_lf = -1 * tr.m3db[0][i - 1] //проводимость слева от узла через нижнюю диагональ
                g_rh = -1 * tr.m3db[2][i + 1] //проводимость справа от узла через верхнюю диагональ
                tr.I[i] = 0.5 * ((tr.U[i - 1] - tr.U[i]) * g_lf + (tr.U[i] - tr.U[i + 1]) * g_rh) //ток  узла как полусумма токов слева и справа от него
                tr.Ignd[i] = tr.Ignd[i - 1] + tr.U[i] * (tr.m3db[1][i] + tr.m3db[0][i - 1] + tr.m3db[2][i + 1]) //ток в земле для узла
            }
            //ток в рельсах для последнего узла
            g_lf = -1 * tr.m3db[2][N - 1] //проводимость слева от последнего узла через верхнюю диагональ последний элемент
            g_rh = 1 / tr.Rvn       //проводимость справа от последнего узла через волновое сопротивление в конце
            tr.I[N - 1] = 0.5 * ((tr.U[N - 2] - tr.U[N - 1]) * g_lf + (tr.U[N - 1] - 0) * g_rh) // ток в рельсе для псоледнего узла
            tr.Ignd[N - 1] = tr.Ignd[N - 2] + tr.U[N - 1] * (tr.m3db[1][N - 1] + tr.m3db[2][N - 1] - 1 / tr.Rvn) //ток в земле для последнего узла
        }
    }
    /**
     * проверка коорректности входного двумерного массива I_poisk[] содержащего токи МПС
     * возвращает true если корректен, и в противном случае false
     */
    private fun verifi_I_poisk(I_poisk: Array<Real>): Boolean {
        if (I_poisk.size != mpss.size) {
            err_and_mes.messeg_solver_error = "Ошибка при вводе массива токов МПС. Количество элементов массива должно совпадать с количеством соответствующих точек МПС"
            return false
        }
        return true
    }

    /**
     *  Возвращает параметр распределенный по узлам сетки в привязке к координатам
     *  @param mesh - сетка пути
     *  @param xRail - массив координат, км;
     *  @param arrParamNode - массив параметра рсапредленный по узлам сетки
     */
    private fun return_X_param_rail(mesh: Mesh, xRail: DoubleArray, arrParamNode: Array<Real> ): Array<Real> {
        //проверки при некорректности возврат null
        if (err_and_mes.data_error || !err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            throw Exception("Расчет не выполнен")
        }
        if (xRail.first() < mesh.X.first() ) { // проверка если координаты расчётного массива за пределами сетки
            throw trackOutOfMeshException("trackOutOfMeshException point ${xRail.first()} out of left edge ${mesh.X.first()} of mesh")
        }
        if (xRail.last() > mesh.X.last()) { // проверка если координаты расчётного массива за пределами сетки
            throw trackOutOfMeshException("trackOutOfMeshException point ${xRail.last()} out of right edge ${mesh.X.last()} of mesh")
        }
        //непосредственно распредление параметра в узлах сетки на массив координат в пределах сетки
        val out = Array<Real>(xRail.size){0.R}
        for (i in xRail.indices) {
            val indexes = mesh.find_2near_index_over_mesh(xRail[i]) // индексы левого и правого узла сетки относительно координаты X_rail[i]
            val proportion_node_left = (mesh.X[indexes[1]] - xRail[i]) / mesh.dX // доля для левого узла сетки обратно пропорциональна расстоянию до правого узла
            out[i] = arrParamNode[indexes[0]] * proportion_node_left + arrParamNode[indexes[1]] * (1 - proportion_node_left) // напряжение в левом узле на его долю + напряжение в правом узле на его долю
        }
        return out
    }
    /**
     * Геттер возвращает напряжение рельс-земля в виде массива по заданным координатам точек
     * @param track - путь
     * @param xRrail - массив значений
     * @return -
     */
    fun get_U_rail( track: Track, xRrail: DoubleArray ): Array<Real> { //  num_track - номер пути; X_rail - массив координат км;
        return return_X_param_rail(track.mesh, xRrail, track.U)
    }
    /**
     * сеттер задает значения в МПС без выполнения непосредственного системы уравнений
     * возвращает true при успешном выполнении процедуры, false  - в противном случае
     * на вход принимает двумерный массив содержащий три строки, это массивы токов в следующем порядке: МПС, ЗАЗ1, ЗАЗ2
     * @param I_poisk - массив токов междупутных соединителей
     */
    fun set_I_poisk_no_call(I_poisk: Array<Real>): Boolean {
        if (err_and_mes.data_error) { // проверка если ошибка исходных данных, то выводим  false. Без кооректного задания ФОТ, ЭПС1 и ЭПС2 - расчёт даже в этом случае невозможен
            err_and_mes.messeg_solver_error = "Задание тока в МПС не возможно. Ошибка исходных данных"
            return false
        }
        if (!verifi_I_poisk(I_poisk)) { //вызываем метод проверки корректности введенного двумерного массива I_poisk
            return false
        }
        //записываем массивы токов в поисковых точках из заданных на входе
        I_mps = I_poisk
        zeros_vector_b() // обнуление вектора правой части всех путей
        eval_U_const() // расчёт от постоянных источников
        eval_node_from_all_I() // расчёт всех узлов с учётом токов МПС
        return true
    }

    /**
     * Расчитывает мгновенную мощность потерь в ОТС, кВт
     * @return мгновенную мощность потерь в ОТС (все элементы ОТС) в кВт
     */
    fun get_P_ots(): Real {
        if (err_and_mes.data_error || !err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим  -1
            return -1.R
            // throw Exception("Расчет не выполнен")
        }
        var P_ots = 0.0.R // мощность ОТС накапливается в этой переменной
        for ( tr in tracks ){
            for ( fot in tr.fot){
                P_ots += -fot.value * get_U_rail(tr, doubleArrayOf(fot.point)).first()//tr.U.first()
            }
            for ( eps in tr.eps){
                P_ots += eps.value * get_U_rail(tr, doubleArrayOf(eps.point)).first()//tr.U.first()
            }
        }
        return P_ots
    }
    /**
     * Проверяет исходные данные расчета
     * @param verify_I условие проверки токов
     * @param verify_Xfot_mps - проверка коорд ФОТ, МПС, ЗАЗ
     * @param verify_Xeps - проверка коорд ЭПС
     */
    private fun verify_data( verify_I: Boolean, verify_Xfot_mps: Boolean, verify_Xeps: Boolean ) {
        err_and_mes.reset_data_error() // очистка data_error
        // проверка по заданию ФОТ и ЭПС и токам
        var Neps_fot_all = 0
        var Isum_fot = 0.0
        var Isum_eps = 0.0
        if (verify_I) {
            // проверяем чтобы сумма токов ЭПС равнялась сумме токов ФОТ не превышая погрешности
            // и что вообще ЭПС  и ФОТ заданы
            for ( tr in tracks ) {
                Neps_fot_all += tr.fot.size + tr.eps.size
                tr.fot.forEach { Isum_fot += it.value.mod }
                tr.eps.forEach { Isum_eps += it.value.mod }
            }
            if (Neps_fot_all == 0) { // если массивы ФОТ и/или ЭПС не заданы
                err_and_mes.data_error = true
                err_and_mes.messeg_data_error += "Массивы ФОТ и/или ЭПС не заданы. "
            }
            val deltaI = Math.abs(Isum_fot - Isum_eps)
            if (deltaI / (Math.max(Math.abs(Isum_fot), Math.abs(Isum_eps)) + 1e-6) > 0.05) {
                err_and_mes.messeg_data_error += " Предупреждение: сумма токов ФОТ и ЭПС расходится более чем на 5 %"
            }
        }

        // проверки по сетке и координатам
        if (verify_Xfot_mps) {
            // проверяем чтобы  границы сетки заданы корректно
            for ( mesh in meshes) {
                if (mesh.X_end - mesh.X_beg <= 2 * mesh.dX) { // если ошибка в границах сетки и менее двух узлов
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "Сетка номер MESH_IDENT границы сетки заданы не корректно, либо получается менее трёх узлов"
                }
            }
            //проверяем чтобы координаты точек ФОТ, ЗАЗ, R_tch укладывались в границы сетки данного пути
            for (tr in tracks){
                for (fot in tr.fot) {
                    if (fot.point > tr.mesh.X_end || fot.point < tr.mesh.X_beg) {
                        err_and_mes.data_error = true
                        err_and_mes.messeg_data_error += "Путь номер ${tr.name}: координаты точек ФОТ выходят за границы сетки"
                    }
                }
                for (zaz in tr.zaz){
                    if (zaz.point > tr.mesh.X_end || zaz.point < tr.mesh.X_beg) {
                        err_and_mes.data_error = true
                        err_and_mes.messeg_data_error += "Путь номер ${tr.name}: координаты точек ЗАЗ выходят за границы сетки"
                    }
                }
                for (tch in tr.Rtch) {
                    if (tch.point > tr.mesh.X_end || tch.point < tr.mesh.X_beg) {
                        err_and_mes.data_error = true
                        err_and_mes.messeg_data_error += "Путь номер ${tr.name}: координаты сосредоченных сопростивлений выходят за границы сетки"
                    }
                }
            }
            // проверка чтобы точки подключения МПС к путям в пределах сетки
            for (mps in mpss) {
                val mesh1 = tracks[mps.startTrack].mesh
                val mesh2 = tracks[mps.endTrack].mesh
                if (mps.startPoint > mesh1.X_end || mps.startPoint < mesh1.X_beg) {
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "МПС ${mps.toString()}: координата начальной точки подключения к пути ${mps.startTrack} выходит за границы сетки"
                }
                if (mps.endPoint > mesh2.X_end || mps.endPoint < mesh2.X_beg) {
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "МПС ${mps.toString()}: координата конечной точки подключения к пути ${mps.endTrack} выходит за границы сетки"
                }
            }
        }
        if (verify_Xeps) {
            //проверяем чтобы координаты точек ЭПС укладывались в границы сетки данного пути
            for ( tr in tracks ) {
                val mesh = tr.mesh
                for ( eps in tr.eps ) {
                    if (eps.point > mesh.X_end || eps.point < mesh.X_beg) {
                        err_and_mes.data_error = true
                        err_and_mes.messeg_data_error += "Путь ${tr.name}: координаты поездов выходят за границы сетки"
                    }
                }
            }
        }
    }
}