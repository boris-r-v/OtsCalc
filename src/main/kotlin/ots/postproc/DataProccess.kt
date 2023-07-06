package ots.postproc

import ots.calc.Compute
import ots.complex.Complex
import ots.complex.mod
import java.util.*
typealias valueDat = ArrayList<ArrayList<Double>>
/**
 *
 */
class TrackStat( track: Int) {
    /*массив номер пути[ массив узлов сетки [ массив значений в этом узле ] ]*/
    internal var U = Array<ArrayList<ArrayList<Double>>>(track) { arrayListOf(arrayListOf()) }
    internal var I = Array<ArrayList<ArrayList<Double>>>(track) { arrayListOf(arrayListOf()) }

    internal var maxU = valueDat()
    internal var minU = valueDat()
    internal var avrU = valueDat()
    internal var skzU = valueDat()
    internal var posAvrU = valueDat()
    internal var rmsI = valueDat()
}
class DataProccess (num_tracks: Int) {
    private var data = TrackStat(num_tracks)
    fun add_one(cmp: Compute) {
        cmp.tracks.forEachIndexed { i, tr ->
            println("U ${tr.name}: ${Arrays.deepToString(tr.U.mod())} ")
            println("I ${tr.name}: ${Arrays.deepToString(tr.I.mod())} ")
            tr.U.forEachIndexed { ii, pnt ->
                println (" Element ${i} ${ii}")
                try {
                    data.U[i].get( ii );
                } catch ( e: IndexOutOfBoundsException ) {
                    println( "Add new U element ${e}" )
                    data.U[i].add( arrayListOf() );
                }
                println ("Add new element ${pnt.mod}")
                data.U[i][ii].add(pnt.mod)
            }
            tr.I.forEachIndexed { ii, pnt ->
                try {
                    data.I[i].get( ii );
                } catch ( e: IndexOutOfBoundsException ) {
                    println( "Add new I element  ${e}" )
                    data.I[i].add( arrayListOf() );
                }
                data.I[i][ii].add(pnt.mod)
            }
        }
    }
    /*массив номер пути[ массив узлов сетки [ массив значений в этом узле ] ]*/
    fun run_post_proccess(){
        data.U.forEachIndexed { i, tr ->
            tr.forEachIndexed { ii, mesh_point ->
                data.maxU[i][ii] = mesh_point.max()
                data.minU[i][ii] = mesh_point.min()
                data.avrU[i][ii] = mesh_point.average()
                data.posAvrU[i][ii] = mesh_point.sumOf { it -> if (it > 0.0) it else 0.0 } / mesh_point.size
                var sd = 0.0
                mesh_point.forEach{ it ->
                    sd += Math.pow(it - data.avrU[i][ii], 2.0 )
                }
                data.skzU[i][ii] = Math.sqrt( sd / mesh_point.size )
            }
        }
        data.I.forEachIndexed{ i, tr ->
            tr.forEachIndexed { ii, mesh_point ->
                data.rmsI[i][ii] = Math.sqrt( mesh_point.sumOf {it*it} / mesh_point.size )
            }
        }
    }
}




