import org.kotlinmath.*
import ots_calc.Mesh

fun main(args: Array<String>) {

    val mesh = Mesh(1.0,5.0, 0.1)
    val x = 3.345

    println("find point: ${x}, " +
            "${mesh.get(mesh.find_near_index_over_mesh(x))}")
    println("find point: ${x}, " +
            "${mesh.find_2near_index_over_mesh(x)[0]} ${mesh.find_2near_index_over_mesh(x)[1]}, " +
            "${mesh.get(mesh.find_2near_index_over_mesh(x))[0]} ${mesh.get(mesh.find_2near_index_over_mesh(x))[1]}")


}