import ots.calc.*
import ots.calc.Mesh
import ots.calc.Track
import ots.complex.*
import ots.statistic.Data
import ots.statistic.MoveImitator

/**
 * -------!!Пример расчета участка ОТС переменного тока с обратным проводом!!----------
 *
 *  Наиболее рационально принимать ОП также как длинную электрическую линию по аналогию с рельсами.
 *  Это означает что для каждого пути с обратным проводом в расчетном классе будет две длинных электрических линии: одна – рельсы, другая – ОП.
 * Никаких функциональных дополнений/изменений к классам модуля расчета ОТС переменного тока Kotlin не требуется.
 * ОТС участка тяговой сети с ОП строится на существующих расчетных классах с существующим функционалом
 * Как и рельсы пути ОП реализуется классом Track. Рельсы и ОП одного пути обязательно должны иметь общую сетку
 *
 */
fun main() {

    /*
    Пример расчета одной мгновенной схемы
   для двухпутного участка 25 кВ с ОП
    */

    /**
     * Создать расчетные сетки.
     * Зададим начальные и конечные координаты и шаг сетки
     */
    val mesh0 = Mesh(8.0,52.0, 0.1)  // общая сетка гл. путей для Рельсов и ОП

    /**
     * Создадим массивы наведенных напряжений вдоль пути в зависимости от координаты
     * Если напряжение меняется то элементов PV может быть любое количество
     */
    val u0 = arrayOf( PV(0.0, 0.0.R) ) // для гл. путей рельсов и ОП - сетка mesh0

    /**
     * Зададим погонное сопротивление рельсов для рельсов и ОП главных путей
     */
    val r_RR = arrayOf( PV(8.0, 0.12615.R+0.5871.I) ) // рельсы
    val r_OP = arrayOf( PV(8.0, 0.1276.R+0.7150.I) ) // ОП
    /**
     * Зададим переходное сопротивление на землю для рельсов и ОП главных путей
     */
    val rp_RR = arrayOf( PV(8.0, 20.0.R) ) // рельсы
    val rp_OP = arrayOf( PV(8.0, 5.0.R) ) // ОП
    /**
     * Зададим расположение и токи фидеров отсоса для данной мгновенной схемы
     */
    val fot0 = arrayOf( PV(10.0, 700.R), PV(50.0, 500.R) ) // отсосы по первому главному пути

    /**
     * Зададим расположение и токи ЭПС
     */
    val eps0 = arrayOf( PV(15.0, 400.R), PV(48.0, 300.R) )  // гл. путь1
    val eps1 = arrayOf( PV(12.0, 300.R), PV(42.0, 350.R))  // гл. путь2


    /**
     * Пустой массив чего-то, для задания отсутствующих элементов в данной мгновенной схеме
     * Необходим технически для отсутствующих массивов в параметрах класса Track
     */
    val emp = arrayOf<PV>()

    /**
     * Создадим пути
     *  имя пути задается исключительно для удобства.
     *  Каждому пути указываем его секту, погонное продольное и переходное сопротивление, расположение отсосов(при наличии) ЭПС (при наличии), наведенные напряжения от КС.
     *  Указываем нужно ли рассчитывать сопротивления по границам сетки или они заданы явно, поля iRv0, iRvn
     */
    val Rv=10000.0.R
    val trackRR1 = Track("RR1", mesh0, r_RR, rp_RR, fot0, eps0, null, null, u0 )      // рельсы гл. путь1
    val trackRR2 = Track("RR2", mesh0, r_RR, rp_RR, emp,  eps1, null, null, u0)      // рельсы гл. путь2
    val trackOP1 = Track("OP1", mesh0, r_OP, rp_OP, emp, emp, null, null, u0 )      // рельсы гл. путь1
    val trackOP2 = Track("OP2", mesh0, r_OP, rp_OP, emp,  emp,  null, null, u0)      // рельсы гл. путь2


    /**
     * Создадим массив междупутных соединителей.
     * Указываем для каждого соединителя:
     * track, которые он соединяет
     * координаты, в которых соединяет track
     * сопротивление соединителя
     */
    var Rmps_RR=0.9e-3.R
    var Rmps_OP=1.5e-3.R
    val mps = arrayOf(
        Mps(trackRR1, trackRR2, 10.0, 10.0, Rmps_RR), /*МПС соединения рельсов гл. путь1 и гл. путь 2*/
        Mps(trackRR1, trackRR2, 30.0, 30.0, Rmps_RR),
        Mps(trackRR1, trackRR2, 50.0, 50.0, Rmps_RR),

        Mps(trackRR1, trackOP1, 10.0, 10.0, Rmps_OP), /*МПС соединения ОП и рельсов гл. пути1*/
        Mps(trackRR1, trackOP1, 15.0, 15.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 20.0, 20.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 25.0, 25.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 30.0, 30.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 35.0, 35.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 40.0, 40.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 45.0, 45.0, Rmps_OP),
        Mps(trackRR1, trackOP1, 50.0, 50.0, Rmps_OP),

        Mps(trackRR2, trackOP2, 10.0, 10.0, Rmps_OP), /*МПС соединения ОП и рельсов гл. пути2*/
        Mps(trackRR2, trackOP2, 15.0, 15.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 20.0, 20.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 25.0, 25.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 30.0, 30.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 35.0, 35.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 40.0, 40.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 45.0, 45.0, Rmps_OP),
        Mps(trackRR2, trackOP2, 50.0, 50.0, Rmps_OP),
    )

    /**
     *    Далее задается междупутное взаимное сопротивление. Оно задается между путями, у которых одна общая сетка.
     *    Т.е. при взаимном сопротивлении фактически должна быть общая система координат вдоль пути.
     *    Поэтому ниже оно задано между всеми тремя главными путями и между двумя путями отходящей двухпутной ветки.
     *    Тупик имеет свою отдельную сетку и не связан с другими путями индуктивно.
     *    Также индуктивно считаются не связанными главные пути и отходящие пути второй ветки
     *    Взаимное сопротивление реализовано через классы RelativeResist и MeshRelativeResist
     */
    /*взаимные сопротивления между разными track*/
    val k_rel=1.0
    val relRR_RR=(0.0493.R+0.3066.I)*k_rel      //между рельсами разных путей
    val relOP_OP=(0.0493.R+0.2466.I)*k_rel     //между ОП разных путей
    val relRR_OP_same=(0.0493.R+0.2918.I)*k_rel      //между рельсами и ОП одного пути
    val relRR_OP_differ=(0.0493.R+0.2615.I)*k_rel      //между рельсами и ОП разных путей
    val relativeResist = RelativeResist()   //объект хранит междупутные сопротивления
    val meshRelativeResist = MeshRelativeResist( mutableMapOf(
        MRRKey(trackRR1, trackRR2) to arrayOf(PV(8.0, relRR_RR)),   // между рельсами путей 1 и 2
        MRRKey(trackOP1, trackOP2) to arrayOf(PV(8.0, relOP_OP)), // между ОП путей 1 и 2
        MRRKey(trackOP1, trackRR1) to arrayOf(PV(8.0, relRR_OP_same)), // между ОП и рельсами путь1
        MRRKey(trackOP2, trackRR2) to arrayOf(PV(8.0, relRR_OP_same)), // между ОП и рельсами путь2
        MRRKey(trackOP1, trackRR2) to arrayOf(PV(8.0, relRR_OP_differ)), // между рельсами путь1 и ОП путь2
        MRRKey(trackOP2, trackRR1) to arrayOf(PV(8.0, relRR_OP_differ))))//он же берется и при сочетании ОП путь2 и Рельсы путь1
    relativeResist.set(meshRelativeResist)

    /**
     * Создадим экземпляр расчетного класса для расчета мгновенной схемы
     * Заданные:
     *      положение и токи электроподвижных единиц,
     *      токи фидеров отсоса
     *      Наведенные в рельсах напряжения от контактной подвески
     * -  задают одну мгновенную схему которую можно рассчитать запустив calcOts()
     */
    println(" --------------Инициализация ОТС. Расчет влиян МПС--------------------------")
    val calc = Compute (
        arrayOf(trackRR1,trackRR2,trackOP1,trackOP2),    /*массив путей*/
        mps,                                             /*массив междупутных соедитнителей*/
        arrayOf(mesh0),                                 /*массив сеток*/
        relativeResist,                                 /*массив межупутных сопротивдления*/
    )

    println(" ---------------------------Расчет мгновенной схемы---------------------------------------")
    println(" ")
    println(" ")
    //calc.computingSettings.isDirectSolver=false  // итерационный решатель для проверки
    calc.calcOts()  // расчет мгновенной схемы
    //вывод на экран модулей напряжений и токов каждого track по сетке
mesh0.X.forEachIndexed{i,value ->
    println(trackRR1.U[i].mod.toString()+", "+trackRR1.I[i].mod.toString()+", "+trackRR2.U[i].mod.toString()+", "+trackRR2.I[i].mod.toString()+", "
    +trackOP1.U[i].mod.toString()+", "+trackOP1.I[i].mod.toString()+", "+trackOP2.U[i].mod.toString()+", "+trackOP2.I[i].mod.toString())
}
/*println((calc.getMpsI()[0]*mps[0].resValue).toString())
println((trackRR1.U[mps[0].startMeshIdx]-trackRR2.U[mps[0].startMeshIdx]).toString())
*/
}