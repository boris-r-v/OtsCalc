package ots_calc

import org.kotlinmath.*
import java.lang.Exception
import kotlin.math.abs
import kotlin.math.max

/**
 * Класс производящий расчет токов и напряжений
 *
 * @param tracks - массив путей
 * @param meshes - массив сеток
 * @param mpss - массив междупутных соединиелей
 * @property mpsCurrent - массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
 * @property numTrackMps - двумерный номеров путей для МПС (в каждой строке начальная и конечная точка)
 * @property indexMps - двумерный массив индексов узлов по сетке МПС (в каждой строке начальная и конечная точка)
 * @property aXFind - матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя
 * @property diffVolt - массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения
 * @property computingSettings
 * @property errorsAndMessages
 */
class Calc(
    private val tracks: Array<Track>,
    private val meshes: Array<Mesh>,
    private val mpss: Array<Mps>,
)
{

    private lateinit var mpsCurrent: Array<Real>
    private lateinit var numTrackMps: Array<IntArray>
    private lateinit var indexMps: Array<IntArray>
    private lateinit var aXFind: Array<Array<Real>>
    private lateinit var diffVolt: Array<Real>
    private val computingSettings: ComputingSettings = ComputingSettings()
    private val errorsAndMessages: ErrorsAndMessages = ErrorsAndMessages()
    init{
        verifyData(verifyI = true, verifyXfotMps = true, verifyXeps = true)
        tracks.forEach { it.m3db = it.mesh.create3diagMatrixBand(it) }
        findIndexMps()
        evalAXFind()
    }

    /**
     * Заполняет массивы свойств
     */
    private fun findIndexMps(){
        numTrackMps = Array(mpss.size) { IntArray(2) } // инициализируем массивы
        indexMps = Array(mpss.size) { IntArray(2) }
        aXFind =Array(mpss.size) { Array(mpss.size){0.R} } // инициализируем общий массив влияния
        for (i in mpss.indices) { //заполняем массивы
            numTrackMps[i][0] = mpss[i].startTrack
            numTrackMps[i][1] = mpss[i].endTrack
            indexMps[i][0] = tracks[mpss[i].startTrack].mesh.findNearIndexOverMesh(mpss[i].startPoint)
            indexMps[i][1] = tracks[mpss[i].endTrack].mesh.findNearIndexOverMesh(mpss[i].endPoint)
        }

    }
    /** Метод для решения СЛАУ с 3диагонал ленточной матрицей
     * методом двойной проходки
     * на входе матрица и вектор правой части
     * на выходе вектор ответов
     */
    private fun solve3diagBand(matrixBand: Array<Array<Real>>, vectorB: Array<Real> ): Array<Real> { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части
        val n = vectorB.size
        val v = Array(n){0.R}
        val u = Array(n){0.R}
        val out = Array(n){0.R}
        //прямая проходка
        v[0] = matrixBand[2][1] / -matrixBand[1][0]
        u[0] = -vectorB[0] / -matrixBand[1][0]
        for (i in 1 until n - 1) {
            v[i] = matrixBand[2][i + 1] / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
            u[i] = (matrixBand[0][i - 1] * u[i - 1] - vectorB[i]) / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
        }
        v[n - 1] = 0.0.R
        u[n - 1] = (matrixBand[0][n - 2] * u[n - 2] - vectorB[n - 1]) / (-matrixBand[1][n - 1] - matrixBand[0][n - 2] * v[n - 2])
        //обратная проходка
        out[n - 1] = u[n - 1]
        for (i in n - 1 downTo 1) {
            out[i - 1] = out[i] * v[i - 1] + u[i - 1]
        }
        return out
    }
    /**
     * Фукция расчета коэффициентов влиияния МПС от них самих
     */
    private fun evalAXFind(){
        val n = this.mpss.size /* количестов поисковых точек */
        val current = 1000.0.R
        // уловный ток МПС для определения коэффициентов влияния
        var u1: Array<Real>
        var u2: Array<Real>
        // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
        val a1XFind = Array(n) { Array(n){0.R} }
        // временные массивы влияния в МПС в начальной и
        val a2XFind = Array(n) { Array(n){0.R} } // в конечной
        // Перебераем поисковые точки МПС
        for (i in this.mpss.indices) {
            tracks[numTrackMps[i][0]].vectorB[indexMps[i][0]] = -current // задаём ток в начальной точке подключения в данном пути
            tracks[numTrackMps[i][1]].vectorB[indexMps[i][1]] = current // задаём ток в конечной точке подключения в данном пути

            u1 = solve3diagBand(
                tracks[numTrackMps[i][0]].m3db,
                tracks[numTrackMps[i][0]].vectorB
            ) //снимаем напряжение по сетке в пути для начальной точки подключения
            u2 = solve3diagBand(
                tracks[numTrackMps[i][1]].m3db,
                tracks[numTrackMps[i][1]].vectorB
            ) //снимаем напряжение по сетке в пути для конечной точки подключения
            tracks[numTrackMps[i][0]].vectorB[indexMps[i][0]] = 0.R // обнуляем ток в текущей начальной точке подключения
            tracks[numTrackMps[i][1]].vectorB[indexMps[i][1]] = 0.R // обнуляем ток в текущей конечной точке подключения
            for (j in this.mpss.indices) { //снова проходим по всем точкам МПС
                if (numTrackMps[i][0] == numTrackMps[j][0]) { // если номер пути начальной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1XFind[i][j] = u1[indexMps[j][0]] / current // находим коэффициенты влияния для текущей точки
                }
                if (numTrackMps[i][1] == numTrackMps[j][0]) { // если номер пути конечной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1XFind[i][j] = u2[indexMps[j][0]] / current // находим коэффициенты влияния для текущей точки
                }
                if (numTrackMps[i][0] == numTrackMps[j][1]) { // если номер пути начальной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2XFind[i][j] = u1[indexMps[j][1]] / current // находим коэффициенты влияния для текущей точки
                }
                if (numTrackMps[i][1] == numTrackMps[j][1]) { // если номер пути конечной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2XFind[i][j] = u2[indexMps[j][1]] / current // находим коэффициенты влияния для текущей точки
                }
                aXFind[i][j] = a1XFind[i][j] - a2XFind[i][j]
            }
        }
    }
    /**
     * Заполняет вектор правой части по заданному двумерному массиву с коэффициентом тока для точек с одним узлом
     * @param mesh - сетка
     * @param vectorB - вектор правой части
     * @param arrXi - массив PV
     * @param coeffI - коэффициент на который умножается значение тока
     */
    private fun valuesVectorB1node(mesh: Mesh, vectorB: Array<Real>, arrXi: Array<PV>, coeffI: Real ): Array<Real>
    {
        for (i in arrXi.indices) {
            vectorB[mesh.findNearIndexOverMesh(arrXi[i].point)] += arrXi[i].value * coeffI
        }
        return vectorB
    }
    /** Добавляет значения в вектор правой части по заданному двумерному массиву с коэффициентом тока для точек с двумя узлами
     * @param mesh - сетка
     * @param vectorB - вектор правой части А
     * @param arrXi -  массив PV
     * @param coeffI  - коэффициент тока
     */
    private fun valuesVectorB2node(mesh: Mesh, vectorB: Array<Real>, arrXi: Array<PV>, coeffI: Real  ): Array<Real>
    {
        for (i in arrXi.indices) {
            val index2 = mesh.find2nearIndexOverMesh( arrXi[i].point ) // номер левого и правого узла
            val k = (arrXi[i].point - mesh.x[0]) / mesh.dX // дробный индекс точки отнсительно номеров узлов сетки
            val current2 = arrXi[i].value * coeffI * (k - index2[0]) // ток в левый узел
            val current1 = arrXi[i].value - current2 // ток в правый узел
            vectorB[index2[0]] += current1
            vectorB[index2[1]] += current2
        }
        return vectorB
    }
    /**
     *  Рассчит напряжение в узлах каждого пути от постоянных источников
     */
    private fun evalDiffVolt() {
        for (tr in tracks){
            tr.vectorB = valuesVectorB1node(tr.mesh, tr.vectorB, tr.fot, (-1.0).R) // в точках ФОТ ток в одном узле с минусом
            tr.vectorB = valuesVectorB2node(tr.mesh, tr.vectorB, tr.eps, 1.0.R ) // в точках ЭПС ток  в двух ближайших узлах
            tr.U = solve3diagBand(tr.m3db, tr.vectorB) //потенциал в рельсах в узлах сетки
            println("tr.name.u, ${tr.U.contentToString()} ")
        }
        diffVolt = Array(mpss.size){ j -> tracks[numTrackMps[j][0]].U[indexMps[j][0]] - tracks[numTrackMps[j][1]].U[indexMps[j][1]]}
        println("U_const:, ${diffVolt.contentDeepToString()} ")
    }

    /**
     * Заполняет вектор правой части всех точек нулями
     */
    private fun zerosVectorB() {
        // проход по всем путям обнуление правой части ФОТ и ЭПС
        for ( tr in tracks ){
            // обнуление для точек в один узел сетки ФОТ
            for (fotPoint in tr.fot ) {
                tr.vectorB[tr.mesh.findNearIndexOverMesh(fotPoint.point)] = 0.0.R
            }
            //обнуление для точек в два узла сетки ЭПС
            for (e in tr.eps) {
                val indx = tr.mesh.find2nearIndexOverMesh(e.point)
                tr.vectorB[indx[0]] = 0.0.R
                tr.vectorB[indx[1]] = 0.0.R
            }
        }
        // проход по всем МПС
        for (i in mpss.indices) {
            tracks[numTrackMps[i][0]].vectorB[indexMps[i][0]] = 0.0.R
            tracks[numTrackMps[i][1]].vectorB[indexMps[i][1]] = 0.0.R
        }
    }
    /**
     * Метод для расчёта мгновенной схемы ОТС с указанием начального значения тока в МПС
     */
    fun calcOts(initIPoisk: Array<Real>): Boolean {
        return if (!verifiIPoisk(initIPoisk)) {
            false
        } else calcIPoisk(initIPoisk)
    }

    /**
     * Метод для расчёта мгновенной схемы ОТС без указания начального значения тока в поисковых точках
     */
    fun calcOts(): Boolean {
        val initIPoisk = Array(mpss.size){0.R}
        return calcIPoisk(initIPoisk)
    }
    /**
     * По сути рассчитывает токи в поисковых точках: МПС. После этого все граничные условия во всех точках с втекающим током определены
     * Расчёт каждого поискового тока сводится к вычислению невязки напряжения на данном элементе и по величине невязки корректируется ток элемента
     * По сути рассчитывается алгебраическая система уравнений методом Ньютона в цикле итераций
     */
    private fun calcIPoisk(initIPoisk: Array<Real>): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        errorsAndMessages.resetSolverError() //обнулим ошибки решателя
        errorsAndMessages.calcCompleted = true
        verifyData(verifyI = false, verifyXfotMps = false, verifyXeps = true) // проверка исходных данных только координ ЭПС, т.к. остальное проверено при инициализации

        if (errorsAndMessages.dataError) { //проверка если данные корректны
            errorsAndMessages.solverError = true
            errorsAndMessages.messegSolverError = "Расчёт невозможен. Ошибка исходных данных"
            return false
        }
        val mpsSize = mpss.size //количество МПС
        val findU = Array(mpsSize){0.R} // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
        val residUMps = Array(mpsSize){0.R}
        var meanResid: Double             // средняя невязка
        val limitMeanResid = computingSettings.convergenceU     // задаём предельную невязку по достижении которой сходимость из класса computing_settings
        var dampingFactor = computingSettings.initialDampingFactor        //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
        var meanResidPred: Double        /* значение невязки по напряжению на предыдущем шаге итераций */
        var iter = 0
        val iterMax = computingSettings.maxIterNumber // счётчик итераций и максимальное число итераций
        var counterNotExceeded: Boolean
        var convergenceNotAchieved: Boolean // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций
        evalDiffVolt() // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counterNotExceeded = true
        convergenceNotAchieved = true
        meanResid = 100.0 //начальное значение средняя невязка до первой итерации
        while (counterNotExceeded && convergenceNotAchieved) {
            meanResidPred = meanResid //предыдущая невязка обновление
            meanResid = 0.0 //текущая невязка скидывается
            for (i in 0 until mpsSize) {
                findU[i] = diffVolt[i] //начинаем с постоянного напряжения (от заданных источников тока ФОТ ЭПС)
                for (j in 0 until mpsSize) {
                    findU[i] += initIPoisk[j] * aXFind[i][j] //добавляем напряжение от МПС
                }
                residUMps[i] = findU[i] - initIPoisk[i] * mpss[i].resValue//невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                initIPoisk[i] += dampingFactor * residUMps[i] / (-0.5 * aXFind[i][i] + mpss[i].resValue) //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования
                meanResid += residUMps[i].mod //обновляем невязку
            }
            meanResid /= mpsSize //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (meanResid > meanResidPred  ) { //FIX ME - тут как то по умному нужно считать невязку в комплексной плоскости
                    dampingFactor *= 0.7
                }
            }
            iter += 1 //обновляем счётчик итераций
            counterNotExceeded = iter < iterMax  // обновляем булевые переменные выхода из цикла итераций
            convergenceNotAchieved = meanResid > limitMeanResid
            //println("iter=$iter mean_resid=$mean_resid damping_factor=$damping_factor")
        }
        computingSettings.currentStateSolver = doubleArrayOf( (iter - 1).toDouble(), meanResid, dampingFactor ) // записываем текущее состояние решателя
        println("I_mps:, ${initIPoisk.contentToString()} ")
        mpsCurrent = initIPoisk // заносим токи в МПС в массивы родительского класса
        evalNodeFromAllI()
        zerosVectorB()
        if (convergenceNotAchieved) {
            errorsAndMessages.solverError = true
            errorsAndMessages.messegSolverError = "Превышено максимальное число итераций равное " + computingSettings.maxIterNumber + " заданная сходимость (сред невязка по напряжен) равная " + computingSettings.convergenceU + " не достигнута."
        }
        return !convergenceNotAchieved // возврат обратное к несходимости: истина если расчёт сошёлся
    }
    /**
     * При найденных токах МПС рассчитывает напряжения и ток в рельсах решением М3ДЛ
     */
    private fun evalNodeFromAllI() {
        var gRh: Real
        var gLf: Real //условная проводимость слева и справа от узла сетки на схеме дискретизации рельсов
        // FIX ME it never usedvar num_mesh: Int  // номер сетки для начального и конечного пути МПС
        var n: Int // количество узлов сетки

        // добавляем в вектор правой часть токи МПС (токи ФОТ и ЭПС добавлены процедурой eval_U_const() )
        for (i in mpss.indices) {
            tracks[numTrackMps[i][0]].vectorB[indexMps[i][0]] -= mpsCurrent[i]
            tracks[numTrackMps[i][1]].vectorB[indexMps[i][1]] += mpsCurrent[i]
        }

        //потенциал в рельсах в узлах сетки для каждого пути
        for (tr in tracks) {
            tr.U = solve3diagBand(tr.m3db, tr.vectorB)
        }

        //расчёт токов в рельсах и тока в земле в узлах сетки для каждого пути
        for (tr in tracks) {
            n = tr.mesh.meshN // число узлов сетки для текущего пути
            //ток в рельсах и земле для первого узла
            gLf = 1 / tr.Rv0 //проводимость слева от первого узла через волновое сопротивление в начале
            gRh = -1 * tr.m3db[0][0] //проводимость справа от первого узла через нижнюю диагональ первый элемент
            tr.I[0] = 0.5 * ((0 - tr.U[0]) * gLf + (tr.U[0] - tr.U[1]) * gRh) //ток первого узла как полусумма токов слева и справа от него
            tr.Ignd[0] = tr.U[0] * (tr.m3db[1][0] + tr.m3db[0][0]) //ток в земле для первого узла
            //ток в рельсах и земле для остальных узлов со второго до предпоследнего
            for (i in 1 until n - 1) {
                gLf = -1 * tr.m3db[0][i - 1] //проводимость слева от узла через нижнюю диагональ
                gRh = -1 * tr.m3db[2][i + 1] //проводимость справа от узла через верхнюю диагональ
                tr.I[i] = 0.5 * ((tr.U[i - 1] - tr.U[i]) * gLf + (tr.U[i] - tr.U[i + 1]) * gRh) //ток  узла как полусумма токов слева и справа от него
                tr.Ignd[i] = tr.Ignd[i - 1] + tr.U[i] * (tr.m3db[1][i] + tr.m3db[0][i - 1] + tr.m3db[2][i + 1]) //ток в земле для узла
            }
            //ток в рельсах для последнего узла
            gLf = -1 * tr.m3db[2][n - 1] //проводимость слева от последнего узла через верхнюю диагональ последний элемент
            gRh = 1 / tr.Rvn       //проводимость справа от последнего узла через волновое сопротивление в конце
            tr.I[n - 1] = 0.5 * ((tr.U[n - 2] - tr.U[n - 1]) * gLf + (tr.U[n - 1] - 0) * gRh) // ток в рельсе для псоледнего узла
            tr.Ignd[n - 1] = tr.Ignd[n - 2] + tr.U[n - 1] * (tr.m3db[1][n - 1] + tr.m3db[2][n - 1] - 1 / tr.Rvn) //ток в земле для последнего узла
        }
    }
    /**
     * проверка коорректности входного двумерного массива I_poisk[] содержащего токи МПС
     * возвращает true если корректен, и в противном случае false
     */
    private fun verifiIPoisk(poiskI: Array<Real>): Boolean {
        if (poiskI.size != mpss.size) {
            errorsAndMessages.messegSolverError = "Ошибка при вводе массива токов МПС. Количество элементов массива должно совпадать с количеством соответствующих точек МПС"
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
    private fun returnXParamRail(mesh: Mesh, xRail: DoubleArray, arrParamNode: Array<Real> ): Array<Real> {
        //проверки при некорректности возврат null
        if (errorsAndMessages.dataError || !errorsAndMessages.calcCompleted) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            throw Exception("Расчет не выполнен")
        }
        if (xRail.first() < mesh.x.first() ) { // проверка если координаты расчётного массива за пределами сетки
            throw TrackOutOfMeshException("trackOutOfMeshException point ${xRail.first()} out of left edge ${mesh.x.first()} of mesh")
        }
        if (xRail.last() > mesh.x.last()) { // проверка если координаты расчётного массива за пределами сетки
            throw TrackOutOfMeshException("trackOutOfMeshException point ${xRail.last()} out of right edge ${mesh.x.last()} of mesh")
        }
        //непосредственно распредление параметра в узлах сетки на массив координат в пределах сетки
        val out = Array(xRail.size){0.R}
        for (i in xRail.indices) {
            val indexes = mesh.find2nearIndexOverMesh(xRail[i]) // индексы левого и правого узла сетки относительно координаты X_rail[i]
            val proportionNodeLeft = (mesh.x[indexes[1]] - xRail[i]) / mesh.dX // доля для левого узла сетки обратно пропорциональна расстоянию до правого узла
            out[i] = arrParamNode[indexes[0]] * proportionNodeLeft + arrParamNode[indexes[1]] * (1 - proportionNodeLeft) // напряжение в левом узле на его долю + напряжение в правом узле на его долю
        }
        return out
    }
    /**
     * Геттер возвращает напряжение рельс-земля в виде массива по заданным координатам точек
     * @param track - путь
     * @param xRrail - массив значений
     * @return -
     */
    private fun getURail(track: Track, xRrail: DoubleArray ): Array<Real> { //  num_track - номер пути; X_rail - массив координат км;
        return returnXParamRail(track.mesh, xRrail, track.U)
    }
    /**
     * сеттер задает значения в МПС без выполнения непосредственного системы уравнений
     * возвращает true при успешном выполнении процедуры, false  - в противном случае
     * на вход принимает двумерный массив содержащий три строки, это массивы токов в следующем порядке: МПС, ЗАЗ1, ЗАЗ2
     * @param poiskI - массив токов междупутных соединителей
     */
    fun setIPoiskNoCall(poiskI: Array<Real>): Boolean {
        if (errorsAndMessages.dataError) { // проверка если ошибка исходных данных, то выводим  false. Без кооректного задания ФОТ, ЭПС1 и ЭПС2 - расчёт даже в этом случае невозможен
            errorsAndMessages.messegSolverError = "Задание тока в МПС не возможно. Ошибка исходных данных"
            return false
        }
        if (!verifiIPoisk(poiskI)) { //вызываем метод проверки корректности введенного двумерного массива I_poisk
            return false
        }
        //записываем массивы токов в поисковых точках из заданных на входе
        mpsCurrent = poiskI
        zerosVectorB() // обнуление вектора правой части всех путей
        evalDiffVolt() // расчёт от постоянных источников
        evalNodeFromAllI() // расчёт всех узлов с учётом токов МПС
        return true
    }

    /**
     * Расчитывает мгновенную мощность потерь в ОТС, кВт
     * @return мгновенную мощность потерь в ОТС (все элементы ОТС) в кВт
     */
    fun getPOts(): Real {
        if (errorsAndMessages.dataError || !errorsAndMessages.calcCompleted) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим  -1
            return (-1).R
            // throw Exception("Расчет не выполнен")
        }
        var otsPower = 0.0.R // мощность ОТС накапливается в этой переменной
        for ( tr in tracks ){
            for ( fot in tr.fot){
                otsPower += -fot.value * getURail(tr, doubleArrayOf(fot.point)).first()//tr.U.first()
            }
            for ( eps in tr.eps){
                otsPower += eps.value * getURail(tr, doubleArrayOf(eps.point)).first()//tr.U.first()
            }
        }
        return otsPower
    }
    /**
     * Проверяет исходные данные расчета
     * @param verifyI условие проверки токов
     * @param verifyXfotMps - проверка коорд ФОТ, МПС, ЗАЗ
     * @param verifyXeps - проверка коорд ЭПС
     */
    private fun verifyData(verifyI: Boolean, verifyXfotMps: Boolean, verifyXeps: Boolean ) {
        errorsAndMessages.resetDataError() // очистка data_error
        // проверка по заданию ФОТ и ЭПС и токам
        var epsFotAll = 0
        var sumFotI = 0.0
        var sumEpsI = 0.0
        if (verifyI) {
            // проверяем чтобы сумма токов ЭПС равнялась сумме токов ФОТ не превышая погрешности
            // и что вообще ЭПС  и ФОТ заданы
            for ( tr in tracks ) {
                epsFotAll += tr.fot.size + tr.eps.size
                tr.fot.forEach { sumFotI += it.value.mod }
                tr.eps.forEach { sumEpsI += it.value.mod }
            }
            if (epsFotAll == 0) { // если массивы ФОТ и/или ЭПС не заданы
                errorsAndMessages.dataError = true
                errorsAndMessages.messegDataError += "Массивы ФОТ и/или ЭПС не заданы. "
            }
            val deltaI = abs(sumFotI - sumEpsI)
            if (deltaI / (max(abs(sumFotI), abs(sumEpsI)) + 1e-6) > 0.05) {
                errorsAndMessages.messegDataError += " Предупреждение: сумма токов ФОТ и ЭПС расходится более чем на 5 %"
            }
        }

        // проверки по сетке и координатам
        if (verifyXfotMps) {
            // проверяем чтобы  границы сетки заданы корректно
            for ( mesh in meshes) {
                if (mesh.endX - mesh.startX <= 2 * mesh.dX) { // если ошибка в границах сетки и менее двух узлов
                    errorsAndMessages.dataError = true
                    errorsAndMessages.messegDataError += "Сетка номер MESH_IDENT границы сетки заданы не корректно, либо получается менее трёх узлов"
                }
            }
            //проверяем чтобы координаты точек ФОТ, ЗАЗ, R_tch укладывались в границы сетки данного пути
            for (tr in tracks){
                for (fot in tr.fot) {
                    if (fot.point > tr.mesh.endX || fot.point < tr.mesh.startX) {
                        errorsAndMessages.dataError = true
                        errorsAndMessages.messegDataError += "Путь номер ${tr.name}: координаты точек ФОТ выходят за границы сетки"
                    }
                }
                for (zaz in tr.zaz){
                    if (zaz.point > tr.mesh.endX || zaz.point < tr.mesh.startX) {
                        errorsAndMessages.dataError = true
                        errorsAndMessages.messegDataError += "Путь номер ${tr.name}: координаты точек ЗАЗ выходят за границы сетки"
                    }
                }
                for (tch in tr.Rtch) {
                    if (tch.point > tr.mesh.endX || tch.point < tr.mesh.startX) {
                        errorsAndMessages.dataError = true
                        errorsAndMessages.messegDataError += "Путь номер ${tr.name}: координаты сосредоченных сопростивлений выходят за границы сетки"
                    }
                }
            }
            // проверка чтобы точки подключения МПС к путям в пределах сетки
            for (mps in mpss) {
                val mesh1 = tracks[mps.startTrack].mesh
                val mesh2 = tracks[mps.endTrack].mesh
                if (mps.startPoint > mesh1.endX || mps.startPoint < mesh1.startX) {
                    errorsAndMessages.dataError = true
                    errorsAndMessages.messegDataError += "МПС ${mps}: координата начальной точки подключения к пути ${mps.startTrack} выходит за границы сетки"
                }
                if (mps.endPoint > mesh2.endX || mps.endPoint < mesh2.startX) {
                    errorsAndMessages.dataError = true
                    errorsAndMessages.messegDataError += "МПС ${mps}: координата конечной точки подключения к пути ${mps.endTrack} выходит за границы сетки"
                }
            }
        }
        if (verifyXeps) {
            //проверяем чтобы координаты точек ЭПС укладывались в границы сетки данного пути
            for ( tr in tracks ) {
                val mesh = tr.mesh
                for ( eps in tr.eps ) {
                    if (eps.point > mesh.endX || eps.point < mesh.startX) {
                        errorsAndMessages.dataError = true
                        errorsAndMessages.messegDataError += "Путь ${tr.name}: координаты поездов выходят за границы сетки"
                    }
                }
            }
        }
    }
}