import org.kotlinmath.*
import ots_calc.Mesh
import ots_calc.PV
import ots_calc.Track
import ots_calc.Calc
import ots_calc.Mps

fun main(args: Array<String>) {

    val mesh0 = Mesh(138.0,180.0, 0.1)
    val mesh1 = Mesh(0.0,7.0, 0.1)
    val mesh2 = Mesh(0.0,14.0, 0.1)

    /*
    val x = 153.345
    println("find point: ${x}, " +
            "${mesh.get(mesh.find_near_index_over_mesh(x))}")
    println("find point: ${x}, " +
            "${mesh.find_2near_index_over_mesh(x)[0]} ${mesh.find_2near_index_over_mesh(x)[1]}, " +
            "${mesh.get(mesh.find_2near_index_over_mesh(x))[0]} ${mesh.get(mesh.find_2near_index_over_mesh(x))[1]}")
*/
    val emp = arrayOf<PV>()
    val r13 = arrayOf( PV(189.0, 0.0254189.R) )
    val rp13 = arrayOf( PV(179.0, 20.R) )

    val r4 = arrayOf( PV(8.0, 0.028.R) )
    val rp4 = arrayOf( PV(8.0, 1.1.R) )

    val r56 = arrayOf( PV(15.0, 0.026.R) )
    val rp56 = arrayOf( PV(15.0, 10.R) )

    val fot0 = arrayOf( PV(140.5, 2300.R), PV(160.2, 2400.R), PV(176.7, 3000.R) )
    val fot4 = arrayOf( PV(10.2, 1400.R) )
    val eps0 = arrayOf( PV(149.0, 800.R), PV(171.2, 3400.R) )
    val eps1 = arrayOf( PV(156.7, 1900.R) )
    val eps2 = arrayOf( PV(145.3, 1600.R) )
    val eps5 = arrayOf( PV(12.1, 1400.R) )

    val track0 = Track("0", r13, rp13, fot0, eps0, emp, emp, mesh0)
    val track1 = Track("1", r13, rp13, emp, eps1, emp, emp,  mesh0)
    val track2 = Track("2", r13, rp13, emp, eps2, emp, emp,  mesh0)
    val track3 = Track("3", r4,  rp4,  emp, emp,  emp, emp,  mesh1)
    val track4 = Track("4", r56, rp56, fot4, emp, emp, emp,  mesh2)
    val track5 = Track("5", r56, rp56, emp, eps5, emp, emp,  mesh2)
    track3.Rv0 = 1e6.R
    track3.Rvn = 1e6.R
    track4.Rv0 = 1e6.R
    track5.Rv0 = 1e6.R

    val mps = arrayOf( Mps(0, 1, 140.5, 140.5, 0.9e-3.R), Mps(1, 2, 140.5, 140.5, 0.9e-3.R),
        Mps(0, 1, 151.5, 151.5, 1.5e-3.R ), /*МПС по главным путям 1-3*/
        Mps(0, 2, 155.5, 155.5, 1.6e-3.R ),
        Mps(0, 1, 160.2, 160.2, 1.4e-3.R ), Mps(1, 2, 160.2, 160.2, 1.1e-3.R ),
        Mps(1, 2, 167.2, 167.2, 1.4e-3.R ),
        Mps(0, 1, 176.7, 176.7, 0.7e-3.R ),
        Mps(1, 2, 177.1, 177.1, 1.8e-3.R ),

        Mps(4, 5, 10.2, 10.2, 1.0e-3.R ), /*МПС по отход2*/
        Mps(0, 3, 152.5, 0.0, 1.0e-5.R ), /*соединение путь гл1 (путь 0 в классе) и однопутн отход тупик (путь 3 в классе)*/
        Mps(1, 4, 170.5, 0.0, 1.0e-5.R ), /*соединение путь гл2 (путь 1 в классе) и  отход2 путь1 (путь 4 в классе) */
        Mps(2, 5, 170.5, 0.0, 1.0e-5.R ), /*соединение путь гл3 (путь 2 в классе) и  отход2 путь2 (путь 5 в классе) */
        )

    val calc = Calc (arrayOf(track0,track1,track2,track3,track4,track5), arrayOf(mesh0,mesh1,mesh2), mps)

    calc.calcOts()

/*
    println("track0.U: ${Arrays.deepToString(track0.U)} ")
    println("track1.U: ${Arrays.deepToString(track1.U)} ")
    println("track2.U: ${Arrays.deepToString(track2.U)} ")
    println("track3.U: ${Arrays.deepToString(track3.U)} ")
    println("track4.U: ${Arrays.deepToString(track4.U)} ")
    println("track5.U: ${Arrays.deepToString(track5.U)} ")

    println("track0.I: ${Arrays.deepToString(track0.I)} ")
    println("track1.I: ${Arrays.deepToString(track1.I)} ")
    println("track2.I: ${Arrays.deepToString(track2.I)} ")
    println("track3.I: ${Arrays.deepToString(track3.I)} ")
    println("track4.I: ${Arrays.deepToString(track4.I)} ")
    println("track5.I: ${Arrays.deepToString(track5.I)} ")

*/
    println("P: ${calc.getPOts()}")

}