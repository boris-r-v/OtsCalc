import org.kotlinmath.*
import ots_calc.Mesh
import ots_calc.PV
import ots_calc.Track
import ots_calc.Calc
import ots_calc.Mps
import java.util.*

fun main(args: Array<String>) {

    val mesh = Mesh(140.0,190.0, 0.1)
/*
    val x = 153.345
    println("find point: ${x}, " +
            "${mesh.get(mesh.find_near_index_over_mesh(x))}")
    println("find point: ${x}, " +
            "${mesh.find_2near_index_over_mesh(x)[0]} ${mesh.find_2near_index_over_mesh(x)[1]}, " +
            "${mesh.get(mesh.find_2near_index_over_mesh(x))[0]} ${mesh.get(mesh.find_2near_index_over_mesh(x))[1]}")
*/

    val r = arrayOf( PV(189.0, 0.0254189.R) )
    val rp = arrayOf( PV(179.0, 20.R) )
    val fot = arrayOf( PV(140.5, 2300.R), PV(160.2, 2400.R), PV(176.7, 3000.R) )
    val eps = arrayOf( PV(149.0, 800.R), PV(171.2, 3400.R) )
    val zaz = arrayOf( PV(179.0, 20.R) )
    val Rtch = arrayOf( PV(179.0, 20.R) )
    val track1 = Track( r,rp,fot,eps,zaz,Rtch, mesh)
    val track2 = Track( r,rp,fot,eps,zaz,Rtch, mesh)

    val mps = arrayOf(Mps(0, 1, 140.5, 140.5, 0.9e-3.R),
                                 Mps(0, 1, 180.5, 180.5, 0.9e-3.R ) )

    val calc = Calc (arrayOf(track1, track2), arrayOf(mesh), mps)

    calc.check()
    calc.calc_ots()


    println("track1.U: ${Arrays.deepToString(track1.U)} ")
    println("track2.U: ${Arrays.deepToString(track2.U)} ")
    println(" ${track1.U.size} ${track2.U.size} " )

    for ( i in 0..track1.U.size-1){
        println("${track1.U[i]} ${track2.U[i]} " )
    }
}