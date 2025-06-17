package ots.calc

import ots.complex.*
import java.lang.Exception
import kotlin.math.abs
import kotlin.math.max

/**
 * Класс производящий расчет токов и напряжений для мгновенной схемы расположения нагрузок
 *
 * @param tracks - массив путей
 * @param mpss - массив междупутных соединиелей
 * @param meshes - массив сеток
 * @param relres - ассоциативный массив взаимных сопротивлений путей принядлежащих одной сетке
 * @property mpsI - массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
 * @property aXFind - матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя
 * @property invaXFind - обратная матрица (к aXFind) коэффициентов влияния тока во всех МПС на напряжения во всех МПС, См.
 * @property knownU - массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения
 * @property computingSettings - настройки расчета
 * @property errorsAndMessages - хранлилище статусов расчета
 */
class Compute(
    internal val tracks: Array<Track>,
    private val mpss: Array<Mps>,
    private val meshes: Array<Mesh>,
    private val relres: RelativeResist,
    )
{
    internal val computingSettings: ComputingSettings = ComputingSettings()
    internal val errorsAndMessages: ErrorsAndMessages = ErrorsAndMessages()

    private val aXFind: Array<Array<Real>> = evalAXFind()
    private val invaXFind: Array<Array<Complex>> = evalInvaXFind()
    private var knownU: Array<Real> = Array(mpss.size) { 0.R }
    private var mpsI: Array<Real> = Array(mpss.size) { 0.R }
    init{
        verifyData(verifyI = true, verifyXfotMps = true, verifyXeps = true) //FIX ME - убрать эту валлидацию в конструкторы объектов с киданием исключений если что-то не так
        /*Установим массив наведенных напряжений из исходных данных после расчета матрицы влияния МПС друг на друга*/
        tracks.forEach { it.setClU() }

//        println("invaXFind=")
//        for (i in 0 until mpss.size) { // предварительно вычтем из матрицы коэффициентов влияния сопротивления МПС
//            println(invaXFind[i].contentDeepToString())
//        }
    }

    /**
     * Метод для решения СЛАУ с 3диагонал ленточной матрицей методом двойной проходки
     * @param matrixBand трехдиагональная матрица
     * @param vectorB вектор свободных членов
     * @param Eks вектор ЭДС по узлам сетки от КС, В/км
     * @param Eots вектор ЭДС по узлам сетки от ОТС, т.е. от других рельсов, В/км
     * @param dX шаг сетки, км
     * @return вектор корнец системы
     */
    private fun solve3diagBand(matrixBand: Array<Array<Real>>, vectorB: Array<Real>, Eks: Array<Real>, Eots: Array<Real>, dX: Double): Array<Real> { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части, Eks и Eots - вектор ЭДС по узлам сетки от КС и ОТС (других рельсов) [В/км], dX- шаг сетки км
        val n = vectorB.size
        val v = Array(n){0.R}
        val u = Array(n){0.R}
        val out = Array(n){0.R}
        var addInode: Complex   //добавка к узловому току - прибавляется к vectorB для данного узла, чтобы учесть ЭДС от КС и ОТС
        //прямая проходка
        addInode=-dX*0.5*(Eots[0]+Eks[0]+Eots[1]+Eks[1])*(-matrixBand[0][0]) //добавка к узловом току 0ой узел =-((1/2)*Eots_0+(1/2)*Eks_0+(1/2)*Eks_1+(1/2)*Eots_1)/r_0   -dX*0.5*0.5*(
        v[0] = matrixBand[2][1] / -matrixBand[1][0]
        u[0] = (-(vectorB[0]+addInode)) / -matrixBand[1][0]
        for (i in 1 until n - 1) {
            addInode=dX*0.5*(Eots[i-1]+Eks[i-1]+Eots[i]+Eks[i])*(-matrixBand[0][i-1]) //первая часть добавки к узловом току iый узел =((1/2)*Eots_i_mn+(1/2)*Eots_i+(1/2)*Eks_i_mn+(1/2)*Eks_i)/r_i_mn   matrixBand[0][i - 1]
            addInode+=-dX*0.5*(Eots[i]+Eks[i]+Eots[i+1]+Eks[i+1])*(-matrixBand[2][i+1]) //вторая часть добавки к узловом току iый узел = -((1/2)*Eots_i+(1/2)*Eks_i+(1/2)*Eks_i_pl+(1/2)*Eots_i_pl)/r_i   matrixBand[2][i + 1]
            v[i] = matrixBand[2][i + 1] / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
            u[i] = (matrixBand[0][i - 1] * u[i - 1] - (vectorB[i]+addInode)) / (-matrixBand[1][i] - matrixBand[0][i - 1] * v[i - 1])
        }
        addInode=dX*0.5*(Eots[n-2]+Eks[n-2]+Eots[n-1]+Eks[n-1])*(-matrixBand[2][n-1]) //добавка к узловом току n-1й узел =((1/2)*Eots_n2+(1/2)*Eots_n1+(1/2)*Eks_n2+(1/2)*Eks_n1)/r_n2    dX*0.5*0.5*(
        v[n - 1] = 0.0.R
        u[n - 1] = (matrixBand[0][n - 2] * u[n - 2] - (vectorB[n - 1]+addInode)) / (-matrixBand[1][n - 1] - matrixBand[0][n - 2] * v[n - 2])
        //обратная проходка
        out[n - 1] = u[n - 1]
        for (i in n - 1 downTo 1) {
            out[i - 1] = out[i] * v[i - 1] + u[i - 1]
        }
        return out
    }

    /**
     * Функция рассчитывает напряжения и токи в узлах группы путей с общей сеткой,
     * предполагается что вектор узловых токов задан и наводимых напряжений от КС тоже.
     * Пересчитываются наводимые от токов в пути этой группы
     * @param mesh сетка в пределах которой считаем наведенные напряжения
     * @param kclU очень важная переменная для расчета, пока пусть будет задел на будущий рефакторинг
     */
    private fun callTracksUnionMesh(mesh: Mesh, kclU: Real=1.0.R ) {
        //FIX ME - kclU -это костыли лучше в будущем заменить
        var avrAbsU = 0.0              // средне напряжение по всем путям группы на текущей итерации
        var avrAbsUold: Double         // тоже самое на итерации назад
        var bIter =true                // признак итерации не исчерпаны
        var bNotConverg =true          // признак сходимость не достигнута
        var k =0                       // счетчик итераций
        var Converg: Double=1.0            // предел относительная сходимость среднего напряжения
        var krlU=computingSettings.initCoefrlU // коэф от наведенного напряжения в track от других track. От 0 до 1. Начинаем с начального в нсатрйоках

        // обнуление всех наводимых напряжений от других путей в группе
        for (tr in mesh.tracks) {
            for ( i in tr.rlU.indices ){
                tr.rlU[i] = 0.0.R
            }
        }
        //расчет напряжений в рельсах и токов при начальных условиях
        for (tr in mesh.tracks ) {
            tr.U = solve3diagBand(tr.m3db, tr.vectorB, tr.clU*kclU, tr.rlU, tr.mesh.dX) // напряжение в рельсах
            evalItrack(tr)    // ток в рельсах
            avrAbsU += tr.U.modAvr()    // добавляем среднее напряжение в данном пути по модулю
        }
        avrAbsUold=avrAbsU                  // принимаем это за значение на прошлой итерации
        // цикл по условию увеличнения коэффициента пропорц от наведенного напряжения k_rlU
        while (krlU.mod<=1.0){
            //println("k_rlU="+k_rlU.mod)
            bNotConverg =true // условие сходимости по среднему абсолютному напряжению
            k = 0 // обнуляем счетчик итерации для вложенного цикла по условию (bIter and bNotConverg)
            // расчет напряжений и токов по условию пока не будет исчерпано число итераций или не достигнута сходимость
            while (bIter and bNotConverg and (k<1)) { // цикл по условию
                // расчет наведенных напряжений вдоль пути от токов в рельсах в группе
                for (tr in mesh.tracks) {
                    //обнуляем у данного пути наводимые напряжения от других путей в группе
                    for (i in tr.rlU.indices) {//обнуляем у данного пути наводимые напряжения от других путей в группе
                        tr.rlU[i] = 0.0.R
                    }
                    //расчитываем у данного пути наводимые напряжения от других путей в группе
                    for (tr2 in mesh.tracks) {
                        if (tr != tr2) {  // условие путь не наводит сам на себя
                            val rlR = relres.get(tr.mesh, tr, tr2) //получение взаимных сопротивлений по сетке
                            tr.rlU = tr.rlU - tr2.I * rlR * krlU // считаем наведенное напряжение от других путей с учетом коэффициента k_rlU
                        }
                    }
                }
                //расчет напряжений в рельсах и токов при текущих условиях
                avrAbsU = 0.0    //среднее значение модулей обнуляем
                for (tr in mesh.tracks) {
                    tr.U = solve3diagBand(tr.m3db, tr.vectorB, tr.clU * kclU, tr.rlU, tr.mesh.dX) // напряжение в рельсах
                    evalItrack(tr)      // ток в рельсах
                    avrAbsU += tr.U.modAvr() // добавляем среднее напряжение в данном пути по модулю
                }
                //заканчиваем цикл обновлением параметров
                Converg =
                    abs(avrAbsU - avrAbsUold) / (avrAbsU + avrAbsUold) // относительная сходимость по среднему абсолютному напряжению на этом и предыдущем шаге
                bNotConverg =
                    (Converg > computingSettings.relativeConvergence)  // признак сходимость не достигнута обновляем
                bIter = (k < computingSettings.maxRelativeIter)                 // признак итерации не исчерпаны  обновляем
                k += 1                                 // счетчик итераций
                avrAbsUold = avrAbsU                  // обновляем значение на прошлой итерации
                }
            // увеличиваем k_rlU на величину в настройках
            krlU = krlU+ computingSettings.stepCoefrlU
        }
        println("k=" + k + ", k_rlU="+(krlU-computingSettings.stepCoefrlU).mod+", avrAbsU=" + avrAbsU + ", Converg=" + Converg)
    }

    /**
     * Функция обнуляет напряжения по и токи всех путей в общей сетке х
     */
    private fun zerosUinMesh(mesh: Mesh) {
        for ( tr in mesh.tracks ){                          //обнуляем массивы напряжений в группе путей пос сетке
            for ( i in tr.U.indices ){
                tr.U[i] = 0.0.R
            }
        }
    }

    /**
     * Функция расчета коэффициентов влияния МПС от них самих
     */
    private fun evalAXFind(): Array<Array<Real>>{
        computingSettings.initCoefrlU=0.009.R
        computingSettings.stepCoefrlU=0.01.R
        computingSettings.relativeConvergence=0.0001
        val out = Array(mpss.size) { Array(mpss.size){0.R} }
        val current = 1000.0.R  // уловный ток МПС для определения коэффициентов влияния
        var u1: Array<Real>
        var u2: Array<Real>
        // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
        val a1XFind = Array(mpss.size) { Array(mpss.size){0.R} }
        // временные массивы влияния в МПС в начальной и
        val a2XFind = Array(mpss.size) { Array(mpss.size){0.R} } // в конечной
        // Перебераем поисковые точки МПС
        mpss.forEachIndexed { i, mps ->
            println()
            println("MPS: 1="+ mps.startTrack.name+", 2="+mps.endTrack.name+", X="+mps.startPoint)
            mps.startTrack.vectorB[mps.startMeshIdx] = -current // задаём ток в начальной точке подключения в данном пути
            mps.endTrack.vectorB[mps.endMeshIdx] = current // задаём ток в конечной точке подключения в данном пути
            callTracksUnionMesh( mps.startTrack.mesh,0.0.R ) //проводим расчёт в группе путей с начальной точкой подключения МПС
            if (mps.startTrack.mesh != mps.endTrack.mesh ){ //если группа путей для конечной точки подключения МПС не совпдает с группой начальной точки
                callTracksUnionMesh( mps.endTrack.mesh,0.0.R ) //проводим расчёт в группе путей с конечной точкой  подключения МПС
            }
            u1 = mps.startTrack.U //снимаем напряжение по сетке в пути для начальной точки подключения МПС
            u2 = mps.endTrack.U  //снимаем напряжение по сетке в пути для конечной точки подключения МПС
            mps.startTrack.vectorB[mps.startMeshIdx] = 0.R // обнуляем ток в текущей начальной точке подключения
            mps.endTrack.vectorB[mps.endMeshIdx] = 0.R // обнуляем ток в текущей конечной точке подключения
            mpss.forEachIndexed { j, mps2 ->
                if (mps.startTrack == mps2.startTrack) { // если номер пути начальной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1XFind[i][j] = u1[mps2.startMeshIdx] / current // находим коэффициенты влияния для текущей точки
                }
                if (mps.endTrack == mps2.startTrack) { // если номер пути конечной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1XFind[i][j] = u2[mps2.startMeshIdx] / current // находим коэффициенты влияния для текущей точки
                }
                if (mps.startTrack == mps2.endTrack) { // если номер пути начальной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2XFind[i][j] = u1[mps2.endMeshIdx] / current // находим коэффициенты влияния для текущей точки
                }
                if (mps.endTrack == mps2.endTrack) { // если номер пути конечной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2XFind[i][j] = u2[mps2.endMeshIdx] / current // находим коэффициенты влияния для текущей точки
                }
                out[i][j] = a1XFind[i][j] - a2XFind[i][j]
            }
            //обнуляем массивы напряжений в группе путей началльной точки подключения МПС
            zerosUinMesh(mps.startTrack.mesh)
            if (mps.startTrack.mesh != mps.endTrack.mesh ){ // если группа путей для конечной точки подключения МПС не совпадает с группой начальной точки
                zerosUinMesh(mps.endTrack.mesh)
            }
        }
        //computingSettings.initCoefrlU=0.049.R
        //computingSettings.stepCoefrlU=0.025.R
        //computingSettings.relativeConvergence=0.001
        return out
    }

    /**
     * Метод просчитывает обратную матрицу к матрице влияний МПС на МПС*
     * @return вектор обратна матрица влияний МПС на МПС
     */
    private fun evalInvaXFind(): Array<Array<Real>>{
        var out = aXFind
        for (i in mpss.indices) { // предварительно вычтем из матрицы коэффициентов влияния сопротивления МПС
            out[i][i] -= mpss[i].resValue
        }
        out = findInvMatGaussJordan(out) // по методу Гаусса-Жордана для комплексных чисел
        return out
    }

    /**
     * Заполняет вектор правой части по заданному двумерному массиву с коэффициентом тока для точек с одним узлом
     * @param mesh - сетка
     * @param vectorB - вектор правой части
     * @param arrXi - массив PV
     * @param coeffI - коэффициент на который умножается значение тока
     */
    private fun valuesVectorB1node(mesh: Mesh, vectorB: Array<Real>, arrXi: Array<PV>, coeffI: Real): Array<Real>
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
    private fun valuesVectorB2node(mesh: Mesh, vectorB: Array<Real>, arrXi: Array<PV>, coeffI: Real): Array<Real>
    {
        for (i in arrXi.indices) {
            val index2 = mesh.find2nearIndexOverMesh( arrXi[i].point ) // номер левого и правого узла
            val k = (arrXi[i].point - mesh.X[0]) / mesh.dX // дробный индекс точки отнсительно номеров узлов сетки
            val current2 = arrXi[i].value * coeffI * (k - index2[0]) // ток в левый узел
            val current1 = arrXi[i].value - current2 // ток в правый узел
            vectorB[index2[0]] += current1
            vectorB[index2[1]] += current2
        }
        return vectorB
    }
    /**
     *  Рассчитывает напряжение в узлах каждого пути от постоянных источников
     */
    private fun evalKnownU() {
        //добавляем в вектор правой части в узлы токи ФОТ и ЭПС по каждому пути
        for (tr in tracks){
            tr.vectorB = valuesVectorB1node(tr.mesh, tr.vectorB, tr.fot, (-1.0).R) // в точках ФОТ ток в одном узле с минусом
            tr.vectorB = valuesVectorB2node(tr.mesh, tr.vectorB, tr.eps, 1.0.R ) // в точках ЭПС ток  в двух ближайших узлах
         }
        // производим расчет напряжений и токов в группах путей по единой сетке
        for (mesh in meshes){
            //println("traks="+mesh.tracks)
            callTracksUnionMesh( mesh)
        }

        //фиксируем напряжения на МПС от заданных источников тока
        knownU = Array(mpss.size){ i -> mpss[i].startTrack.U[mpss[i].startMeshIdx] - mpss[i].endTrack.U[mpss[i].endMeshIdx] }

        println("U_const:, ${knownU.contentDeepToString()} ")
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
        for (mps in mpss) {
            mps.startTrack.vectorB[mps.startMeshIdx] = 0.0.R
            mps.endTrack.vectorB[mps.endMeshIdx] = 0.0.R
        }
    }
    /**
     * Метод для расчёта мгновенной схемы ОТС с указанием начального значения тока в МПС
     */
    fun calcOts(initMpsI: Array<Real>): Boolean {
        return if (!verifiIPoisk(initMpsI)) {
            false
        } else calcIPoisk(initMpsI)
    }

    /**
     * Метод для расчёта мгновенной схемы ОТС без указания начального значения тока в поисковых точках
     * После расчета производится копирование текущих значений U, I в узлах сетки в массивы исторических значений
     */
    fun calcOts(): Boolean {
        var ret: Boolean
        if (computingSettings.isDirectSolver) {  // прямой решатель токов МПС
            ret = calcIPoisk_invMat()
        }
        else {   // итерационный решатель токов МПС
            ret = calcIPoisk(Array(mpss.size){0.R})
        }
        tracks.forEach { it.copy2Hist() }
        return ret
    }

    private fun сheckBeforeIPoisk() : Boolean {
        errorsAndMessages.resetSolverError() //обнулим ошибки решателя
        errorsAndMessages.calcCompleted = true
        verifyData(verifyI = false, verifyXfotMps = false, verifyXeps = true) // проверка исходных данных только координ ЭПС, т.к. остальное проверено при инициализации

        if (errorsAndMessages.dataError) { //проверка если данные корректны
            errorsAndMessages.solverError = true
            errorsAndMessages.messegSolverError = "Расчёт невозможен. Ошибка исходных данных"
            return false
        }
        return true
    }
    /**
     * По сути рассчитывает токи в поисковых точках: МПС. После этого все граничные условия во всех точках с втекающим током определены
     * Расчёт каждого поискового тока сводится к вычислению невязки напряжения на данном элементе и по величине невязки корректируется ток элемента
     * По сути рассчитывается алгебраическая система уравнений методом Ньютона в цикле итераций
     */
    private fun calcIPoisk(initIPoisk: Array<Real>): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        if (!сheckBeforeIPoisk()) {
            return false
        }
        val findU = Array(mpss.size){0.R} // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
        val residUMps = Array(mpss.size){0.R}
        var meanResid: Double             // средняя невязка
        val limitMeanResid = computingSettings.convergenceU     // задаём предельную невязку по достижении которой сходимость из класса computing_settings
        var dampingFactor = computingSettings.initialDampingFactor        //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
        var meanResidPred: Double        /* значение невязки по напряжению на предыдущем шаге итераций */
        var iter = 0
        val iterMax = computingSettings.maxIterNumber // счётчик итераций и максимальное число итераций
        var counterNotExceeded: Boolean
        var convergenceNotAchieved: Boolean // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций
        evalKnownU() // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counterNotExceeded = true
        convergenceNotAchieved = true
        meanResid = 100.0 //начальное значение средняя невязка до первой итерации
        while (counterNotExceeded && convergenceNotAchieved) {
            meanResidPred = meanResid //предыдущая невязка обновление
            meanResid = 0.0 //текущая невязка скидывается
            for (i in mpss.indices ) {
                findU[i] = knownU[i] //начинаем с постоянного напряжения (от заданных источников тока ФОТ ЭПС)
                for (j in mpss.indices) {
                    findU[i] += initIPoisk[j] * aXFind[i][j] //добавляем напряжение от МПС
                }
                residUMps[i] = findU[i] - initIPoisk[i] * mpss[i].resValue//невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                initIPoisk[i] += dampingFactor * residUMps[i] / (-0.5 * aXFind[i][i] + mpss[i].resValue) //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования
                meanResid += residUMps[i].mod //обновляем невязку
            }
            meanResid /= mpss.size //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (meanResid > meanResidPred  ) { //FIX ME - тут как то по умному нужно считать невязку в комплексной плоскости
                    dampingFactor *= 0.7
                }
            }
            iter += 1 //обновляем счётчик итераций
            counterNotExceeded = iter < iterMax  // обновляем булевые переменные выхода из цикла итераций
            convergenceNotAchieved = meanResid > limitMeanResid
            //println("iter=$iter mean_resid=$meanResid damping_factor=$dampingFactor")
        }
        computingSettings.currentStateSolver = doubleArrayOf( (iter - 1).toDouble(), meanResid, dampingFactor ) // записываем текущее состояние решателя
        mpsI = initIPoisk // заносим токи в МПС в массивы родительского класса
        println("I_mps:, ${mpsI.contentToString()} ")
        evalNodeFromAllI()
        if (convergenceNotAchieved) {
            errorsAndMessages.solverError = true
            errorsAndMessages.messegSolverError = "Превышено максимальное число итераций равное " + computingSettings.maxIterNumber + " заданная сходимость (сред невязка по напряжен) равная " + computingSettings.convergenceU + " не достигнута."
        }
        //clarificationIPoisk()
        zerosVectorB()
        return !convergenceNotAchieved // возврат обратное к несходимости: истина если расчёт сошёлся
    }

    /**
     * По сути рассчитывает токи в поисковых точках: МПС через обратную матрицу. По сути прямой решатель
     */
    private fun calcIPoisk_invMat(): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        if (!сheckBeforeIPoisk()) {
            return false
        }
        evalKnownU() // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        for (i in mpss.indices) {
            mpsI[i] = 0.R // предварительно ток МПС =0, чтобы начать с нуля
            for (j in mpss.indices) {  //умножение строки на столбец
                mpsI[i] -= invaXFind[i][j] * knownU[j] //  this.I_mps[i] += -this.inv_a_x_find[j][i]*this.U_const[j];
                }
        }
        println("I_mps:, ${mpsI.contentToString()} ")
        evalNodeFromAllI()
        zerosVectorB()
        return true // возврат обратное к несходимости: истина если расчёт сошёлся
    }

    /**
     * Уточнение IPoisk     *
     */
    private fun clarificationIPoisk(): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        val findU = Array(mpss.size){0.R} // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
        val residUMps = Array(mpss.size){0.R}
        var meanResid: Double             // средняя невязка
        val limitMeanResid = 0.01 //computingSettings.convergenceU     // задаём предельную невязку по достижении которой сходимость из класса computing_settings
        var dampingFactor = 0.1       //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
        var meanResidPred: Double        /* значение невязки по напряжению на предыдущем шаге итераций */
        var iter = 0
        val iterMax = 5 // счётчик итераций и максимальное число итераций
        var counterNotExceeded: Boolean
        var convergenceNotAchieved: Boolean // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций
        val oldInitCoefrlU=computingSettings.initCoefrlU

        println("-------------Уточнение поисковых токов-------------------")
        println("0  mpsI="+mpsI.mod().contentDeepToString())
        //computingSettings.initCoefrlU=0.999.R
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counterNotExceeded = true
        convergenceNotAchieved = true
        meanResid = 100.0 //начальное значение средняя невязка до первой итерации
        while (counterNotExceeded && convergenceNotAchieved) {
            meanResidPred = meanResid //предыдущая невязка обновление
            meanResid = 0.0 //текущая невязка скидывается
            for (i in mpss.indices ) {
                findU[i] = -mpss[i].endTrack.U[mpss[i].endMeshIdx]+mpss[i].startTrack.U[mpss[i].startMeshIdx] //
                residUMps[i] = findU[i] - mpsI[i] * mpss[i].resValue//невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                mpsI[i] += dampingFactor * residUMps[i] / (-10 * aXFind[i][i] + mpss[i].resValue) //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования
                meanResid += residUMps[i].mod //обновляем невязку
            }
            println(iter.toString()+"  residUMps="+residUMps.mod().contentDeepToString())
            meanResid /= mpss.size //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (meanResid > meanResidPred  ) { //FIX ME - тут как то по умному нужно считать невязку в комплексной плоскости
                    dampingFactor *= 0.7
                }
            }
            iter += 1 //обновляем счётчик итераций
            counterNotExceeded = iter < iterMax  // обновляем булевые переменные выхода из цикла итераций
            convergenceNotAchieved = meanResid > limitMeanResid
            println(iter.toString()+"  mpsI="+mpsI.mod().contentDeepToString())
            evalNodeFromAllI() // просчитываем токи
            //println("iter=$iter mean_resid=$mean_resid damping_factor=$damping_factor")
        }
        computingSettings.currentStateSolver = doubleArrayOf( (iter - 1).toDouble(), meanResid, dampingFactor ) // записываем текущее состояние решателя
        //println("I_mps:, ${initIPoisk.contentToString()} ")
        if (convergenceNotAchieved) {
            errorsAndMessages.solverError = true
            errorsAndMessages.messegSolverError = "Превышено максимальное число итераций равное " + computingSettings.maxIterNumber + " заданная сходимость (сред невязка по напряжен) равная " + computingSettings.convergenceU + " не достигнута."
        }
        computingSettings.initCoefrlU=oldInitCoefrlU
        return !convergenceNotAchieved // возврат обратное к несходимости: истина если расчёт сошёлся
    }
    /**
     * При найденных токах МПС рассчитывает напряжения и ток в рельсах решением М3ДЛ
     */
    private fun evalNodeFromAllI() {
         // добавляем в вектор правой часть токи МПС (токи ФОТ и ЭПС добавлены процедурой eval_U_const() )
        mpss.forEachIndexed { i, mps ->
                mps.startTrack.vectorB[mps.startMeshIdx] -= mpsI[i]
                mps.endTrack.vectorB[mps.endMeshIdx] += mpsI[i]
        }
        //потенциал и токи в рельсах в рельсах для всех групп путей с общими сетками
        meshes.forEach { callTracksUnionMesh( it ) }
         //расчёт токов в земле в узлах сетки для каждого пути
        for (tr in tracks) {
            evalIgndTrack(tr) // ток в земле
        }
    }

    /**
     * Рассчитывает ток в узлах сетки, предполагается что напряжения в узлах уже определенны
     * @param track - массив путей
     */
    private fun evalItrack(track: Track) {
        val n = track.mesh.size // число узлов сетки для текущего пути
        var gRh: Real //условная проводимость слева и справа от узла сетки на схеме дискретизации рельсов, См
        var gLf: Real
        var IRh: Real //условный ток слева и справа от узла сетки на схеме дискретизации рельсов, А
        var ILf: Real


        //ток в рельсах для первого узла
        gLf = 1 / track.Rv0 //проводимость слева от первого узла через волновое сопротивление в начале
        gRh = -1 * track.m3db[0][0] //проводимость справа от первого узла через нижнюю диагональ первый элемент
        ILf=(0 - track.U[0]) * gLf // ток слева от узла А
        IRh=(track.U[0] - track.U[1]+(0.5*track.mesh.dX)*(track.clU[0]+track.clU[1]+track.rlU[0]+track.rlU[1])) * gRh // ток справа от узла
        track.I[0] = 0.5*(ILf  + IRh) //ток первого узла как полусумма токов слева и справа от него
         //ток в рельсах для остальных узлов со второго до предпоследнего
        for (i in 1 until n - 1) {
            gLf = -1 * track.m3db[0][i-1] //проводимость слева от узла через нижнюю диагональ
            gRh = -1 * track.m3db[2][i+1] //проводимость справа от узла через верхнюю диагональ
            ILf=(track.U[i - 1] - track.U[i]+0.5*track.mesh.dX*(track.clU[i-1]+track.clU[i]+track.rlU[i-1]+track.rlU[i])) * gLf // ток слева от узла А (u_i_mn-u_i+(1/2*(Eks_i_mn+Eks_i))*dx+(1/2*(Eots_i_mn+Eots_i))*dx)/r_i_mn
            IRh=(track.U[i] - track.U[i+1]+0.5*track.mesh.dX*(track.clU[i]+track.clU[i+1]+track.rlU[i]+track.rlU[i+1])) * gRh // ток справа от узла  (u_i-u_i_pl+(1/2*(Eks_i+Eks_i_pl))*dx+(1/2*(Eots_i+Eots_i_pl))*dx)/r_i
            track.I[i] =0.5*(ILf  + IRh) //ток  узла как полусумма токов слева и справа от него
           }
        //ток в рельсах для последнего узла
        gLf = -1 * track.m3db[2][n - 1] //проводимость слева от последнего узла через верхнюю диагональ последний элемент
        gRh = 1 / track.Rvn       //проводимость справа от последнего узла через волновое сопротивление в конце
        ILf=(track.U[n - 2] - track.U[n-1]+(0.5*track.mesh.dX)*(track.clU[n-2]+track.clU[n-1]+track.rlU[n-2]+track.rlU[n-1])) * gLf // ток слева от узла А
        IRh=(track.U[n-1] - 0) * gRh // ток справа от узла
        track.I[n-1] =0.5*(ILf  + IRh) //ток  узла как полусумма токов слева и справа от него
        }

    /**
     * Рассчитывает ток в земле в узлах сетки, предполагается что напряжения в узлах уже определенны
     * @param track - путь
     */
    private fun evalIgndTrack(track: Track) {
        val n = track.mesh.size // число узлов сетки для текущего пути
        //ток в земле для первого узла
        track.Ignd[0] = track.U[0] * (track.m3db[1][0] + track.m3db[0][0])
        //ток в земле для остальных узлов со второго до предпоследнего
        for (i in 1 until n - 1) {
            track.Ignd[i] = track.Ignd[i - 1] + track.U[i] * (track.m3db[1][i] + track.m3db[0][i - 1] + track.m3db[2][i + 1])
        }
        //ток в земле для последнего узла
        track.Ignd[n - 1] = track.Ignd[n - 2] + track.U[n - 1] * (track.m3db[1][n - 1] + track.m3db[2][n - 1] - 1 / track.Rvn)
    }

    /**
     * Проверка корректности входного двумерного массива I_poisk[] содержащего токи МПС
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
     *  @param arrParamNode - массив параметра распределенный по узлам сетки
     */
    private fun returnXParamRail(mesh: Mesh, xRail: DoubleArray, arrParamNode: Array<Real> ): Array<Real> {
        //проверки при некорректности возврат null
        if (errorsAndMessages.dataError || !errorsAndMessages.calcCompleted) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            throw Exception("Расчет не выполнен")
        }
        if (xRail.first() < mesh.X.first() ) { // проверка если координаты расчётного массива за пределами сетки
            throw TrackOutOfMeshException("trackOutOfMeshException point ${xRail.first()} out of left edge ${mesh.X.first()} of mesh")
        }
        if (xRail.last() > mesh.X.last()) { // проверка если координаты расчётного массива за пределами сетки
            throw TrackOutOfMeshException("trackOutOfMeshException point ${xRail.last()} out of right edge ${mesh.X.last()} of mesh")
        }
        //непосредственно распредление параметра в узлах сетки на массив координат в пределах сетки
        val out = Array(xRail.size){0.R}
        for (i in xRail.indices) {
            val indexes = mesh.find2nearIndexOverMesh(xRail[i]) // индексы левого и правого узла сетки относительно координаты X_rail[i]
            val proportionNodeLeft = (mesh.X[indexes[1]] - xRail[i]) / mesh.dX // доля для левого узла сетки обратно пропорциональна расстоянию до правого узла
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
        mpsI = poiskI
        zerosVectorB() // обнуление вектора правой части всех путей
        evalKnownU() // расчёт от постоянных источников
        evalNodeFromAllI() // расчёт всех узлов с учётом токов МПС
        return true
    }

    /**
     * Рассчитывает мгновенную мощность потерь в ОТС, кВт
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
                val mesh1 = mps.startTrack.mesh
                val mesh2 = mps.endTrack.mesh
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

    /**
     * Возвращает массив токов междупутных соединителей для данной мгновенной схемы
     *@return Массив токов междупутных соединителей
     */
    fun getMpsI(): Array<Real> {
        return mpsI
    }


/*
    //--------------------------------функции для комплексеных массивов (костыли)-------------------------------------------
    /**
     * Сложение комплексных массивов
     *  @param arr1 комплексеный массив1
     *  @param arr2 комплексеный массив2
     * @return Массив суммарный поэлементный
     */
    private fun sumComplexArray(arr1: Array<Real>,arr2: Array<Real> ): Array<Real>{
        val out: Array<Real>
        val n : Int = arr1.size
        out=Array(n) { 0.R }
        for (i in 0 until n - 1) {
            out[i]=arr1[i]+arr2[i]
        }
        return out
    }
    /**
     * Сложение комплексного массива и числа одного
     *  @param arr1 комплексеный массив1
     *  @param OneValue комплексеное число одно
     * @return Массив суммарный поэлементный
     */
    private fun sumComplexArray(arr1: Array<Real>, OneValue: Real ): Array<Real>{
        val out: Array<Real>
        val n : Int = arr1.size
        out=Array(n) { 0.R }
        for (i in 0 until n - 1) {
            out[i]=arr1[i]+OneValue
        }
        return out
    }
    /**
     * умножениее комплексных массивов
     * @param arr1 комплексеный массив1
     * @param arr2 комплексеный массив2
     * @return Массив произведений поэлементный
     */
    private fun mulComplexArray(arr1: Array<Real>,arr2: Array<Real> ): Array<Real>{
        val out: Array<Real>
        val n : Int = arr1.size
        out=Array(n) { 0.R }
        for (i in 0 until n - 1) {
            out[i]=arr1[i]*arr2[i]
        }
        return out
    }
    /**
     * Умножение комплексного массива и числа одного
     * @param arr1 комплексеный массив1
     * @param OneValue комплексеное число одно
     * @return Массив суммарный поэлементный
     */
    private fun mulComplexArray(arr1: Array<Real>, OneValue: Real ): Array<Real>{
        val out: Array<Real>
        val n : Int = arr1.size
        out=Array(n) { 0.R }
        for (i in 0 until n - 1) {
            out[i]=arr1[i]*OneValue
        }
        return out
    }

    /**
     * Суммированние всех элементов массива по модулю
     * @param arr комплексный массив
     * @return сумма модулей элементов массива
     */
    private fun sumAbsElementComlexArray(arr: Array<Real>): Double{
        var out: Double
        val n : Int = arr.size
        out=0.0
        for (i in 0 until n - 1) {
            out=out+arr[i].mod
        }
        return out
    }
*/
}


