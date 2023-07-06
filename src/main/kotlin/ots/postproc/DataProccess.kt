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
    fun handle_one(cmp: Compute) {
        cmp.tracks.forEachIndexed { i, tr ->
            println("U ${tr.name}: ${Arrays.deepToString(tr.U.mod())} ")
            println("I ${tr.name}: ${Arrays.deepToString(tr.I.mod())} ")
            tr.U.forEachIndexed { ii, pnt ->
                data.U[i][ii].add(pnt.mod)
                if (data.maxU[i][ii] < pnt.mod) {
                    data.maxU[i][ii] = pnt.mod
                }
                if (data.minU[i][ii] > pnt.mod) {
                    data.minU[i][ii] = pnt.mod
                }
                data.avrU[i][ii] = data.U[i][ii].sumOf { it } / data.U[i][ii].size
                data.posAvrU[i][ii] = data.U[i][ii].sumOf{ it -> if (it > 0) it else 0 } / data.U[i][ii].size
            }
            tr.I.forEachIndexed { ii, pnt ->
                data.I[i][ii].add(pnt.mod)
            }
        }
    }
}




