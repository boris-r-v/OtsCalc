import ots.calc.*
import ots.calc.Mesh
import ots.calc.Track
import ots.complex.*
import java.util.*

/*
Пример расчета одной мгновенной схемы без описания особого. Для тестирования старый вариант
 */

fun main(args: Array<String>) {

    val mesh0 = Mesh(138.0,180.0, 0.1)
    val mesh1 = Mesh(0.0,7.0, 0.1)
    val mesh2 = Mesh(0.0,14.0, 0.1)

    val u0 = arrayOf( PV(0.0, 0.0.R) )
    val u1 = arrayOf( PV(0.0, 0.0.R ))
    val u2 = arrayOf( PV(0.0, 0.0.R) )

    val r13 = arrayOf( PV(189.0, 0.12615.R+0.5871.I) )
    val rp13 = arrayOf( PV(179.0, 20.R) )

    val r4 = arrayOf( PV(8.0, 0.12615.R+0.5871.I) )
    val rp4 = arrayOf( PV(8.0, 1.1.R) )

    val r56 = arrayOf( PV(15.0, 0.12615.R+0.5871.I) )
    val rp56 = arrayOf( PV(15.0, 10.R) )

    val fot0 = arrayOf( PV(140.5, 2300.R), PV(160.2, 2400.R), PV(176.7, 3000.R) )
    val fot4 = arrayOf( PV(10.2, 1400.R) )
    val eps0 = arrayOf( PV(149.0, 800.R), PV(171.2, 3400.R) )
    val eps1 = arrayOf( PV(156.7, 1900.R) )
    val eps2 = arrayOf( PV(145.3, 1600.R) )
    val eps5 = arrayOf( PV(12.1, 1400.R) )

    val emp = arrayOf<PV>()
    /**
     *     Если поставить сопротивление 0.0.R то расчет совпадает с постоянным током если пос
     *     Если поставить 1мкОм - то тоже почти совпадает
     */
    val track0 = Track("0", mesh0, r13, rp13, fot0, eps0, null, null, u0 )
    val track1 = Track("1", mesh0, r13, rp13, emp,  eps1,  null, null, u0)
    val track2 = Track("2", mesh0, r13, rp13, emp,  eps2, null, null, u0)
    val track3 = Track("3", mesh1, r4,  rp4,  emp,  emp,  1e6.R, 1e6.R, u1 )
    val track4 = Track("4", mesh2, r56, rp56, fot4, emp,  1e6.R, null, u2 )
    val track5 = Track("5", mesh2, r56, rp56, emp,  eps5, 1e6.R, null, u2 )

    val mps = arrayOf(
        Mps(track0, track1, 140.5, 140.5, 0.9e-3.R), /*МПС по главным путям 1-3*/
        Mps(track1, track2, 140.5, 140.5, 0.9e-3.R),
        Mps(track0, track1, 151.5, 151.5, 1.5e-3.R ),
        Mps(track0, track2, 155.5, 155.5, 1.6e-3.R ),
        Mps(track0, track1, 160.2, 160.2, 1.4e-3.R ),
        Mps(track1, track2, 160.2, 160.2, 1.1e-3.R ),
        Mps(track1, track2, 167.2, 167.2, 1.4e-3.R ),
        Mps(track0, track1, 176.7, 176.7, 0.7e-3.R ),
        Mps(track1, track2, 177.1, 177.1, 1.8e-3.R ),

        Mps(track4, track5, 10.2, 10.2, 1.0e-3.R ), /*МПС по отход2*/
        Mps(track0, track3, 152.5, 0.0, 1.0e-5.R ), /*соединение путь гл1 (путь 0 в классе) и однопутн отход тупик (путь 3 в классе)*/
        Mps(track1, track4, 170.5, 0.0, 1.0e-5.R ), /*соединение путь гл2 (путь 1 в классе) и  отход2 путь1 (путь 4 в классе) */
        Mps(track2, track5, 170.5, 0.0, 1.0e-5.R ), /*соединение путь гл3 (путь 2 в классе) и  отход2 путь2 (путь 5 в классе) */
        )
    val rpRes=0.04.R+0.3.I      //оставил тоже значение междупутного соединения
    val rr = RelativeResist()   //объект хранит междупутные сопротивления
    val trr = MeshRelativeResist( mutableMapOf(   MRRKey(track0, track1) to arrayOf(PV(138.0, rpRes)),
                                                            MRRKey(track0, track2) to arrayOf(PV(138.0, rpRes)),
                                                            MRRKey(track1, track2) to arrayOf(PV(138.0, rpRes))
                                                        ))
    rr.set(trr)
    rr.set(MeshRelativeResist( mutableMapOf( MRRKey(track4, track5) to arrayOf(PV(0.0, rpRes)))))
    val calc = Compute (arrayOf(track0,track1,track2,track3,track4,track5),     /*массив путей*/
                        mps,                                                    /*массив междупутных соедитнителей*/
                        arrayOf(mesh0,mesh1, mesh2),                            /*массив сеток*/
                        rr,                                                     /*массив межупутных сопротивдления*/
                        )  // добавил в аргументы массив сеток
    calc.calcOts()
    println(calc.computingSettings.currentStateSolver[0])
    println(calc.computingSettings.currentStateSolver[1])
    println(calc.computingSettings.currentStateSolver[2])
    println("track0.I: ${Arrays.deepToString(track0.I.mod())} ")
    println("track1.I: ${Arrays.deepToString(track1.I.mod())} ")
    println("track2.I: ${Arrays.deepToString(track2.I.mod())} ")
    println(Arrays.deepToString(mesh0.X))
//    println("track3.U: ${Arrays.deepToString(track3.U)} ")
//    println("track4.U: ${Arrays.deepToString(track4.U)} ")
//    println("track5.U: ${Arrays.deepToString(track5.U)} ")
//
//    println("track0.I: ${Arrays.deepToString(track0.I)} ")
//    println("track1.I: ${Arrays.deepToString(track1.I)} ")
//    println("track2.I: ${Arrays.deepToString(track2.I)} ")
//    println("track3.I: ${Arrays.deepToString(track3.I)} ")
//    println("track4.I: ${Arrays.deepToString(track4.I)} ")
//    println("track5.I: ${Arrays.deepToString(track5.I)} ")

}