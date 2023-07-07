package ots.postproc

import ots.calc.Track
import java.util.*

/**
 * Структура данных для хранения статистичесих данных путей в узлах сетки
 * Например максимальные значения U/I в каждом узле сетки
 * Стуктура массив_путей[ массив узлов сетки [ Значение ] ]
 */
typealias meshPointData = Array<ArrayList<Double>>
/**
 * Статистика по одному пути, данные распределены по узлам сетки пути
 * @param mesh_size колическво узлов сетки
 * @param rmsWindow количество шагов расчета за который расчитывать RMS
 * @property maxU массив максимальных значений напряжений в узлах сетки
 * @property minU массив минимальны значений напряжений в узлах сетки
 * @property avrU массив средних арифметических значений напряжений в узлах сетки
 * @property skzU массив среднеквадратичных отклонений от среднего значения напряжений у узлах сетки
 * @property posAvrU массив положительных средних арифметических значений напряжений в узлах сетки
 * @property rmsI действующее значение тока за rms_steps шагов расчета в каждом узле сетки
 */
class TrackStat(track: Track, rmsWindow: Int) {
    internal var maxU: Array<Double> = Array(track.mesh.size) {0.0}
    internal var minU: Array<Double> = Array(track.mesh.size) {0.0}
    internal var avrU: Array<Double> = Array(track.mesh.size) {0.0}
    internal var skoU: Array<Double> = Array(track.mesh.size) {0.0}
    internal var posAvrU: Array<Double> = Array(track.mesh.size) {0.0}
    internal var rmsI: Array<Double> = Array(track.mesh.size) {0.0}

    /**
     * В методе инициализации проходит вся магия приготовления/расчета статистики
     */
    init {
        track.histU.forEachIndexed { i, arrForOneMeshPoint ->
            var sumOfUArray = 0.0
            var sumOfPosUArray = 0.0
            var min = 10000.0
            var max =-10000.0
            var size = arrForOneMeshPoint.size
            arrForOneMeshPoint.forEach{ real ->
                val rmod = real.mod
                if ( min > rmod )
                    min = rmod
                if ( max < rmod )
                    max = rmod
                sumOfUArray += rmod
                if (rmod > 0)
                    sumOfPosUArray += rmod

            }
            maxU[i] = max
            minU[i] = min
            avrU[i] = sumOfUArray / size
            posAvrU[i] = sumOfPosUArray / size
            val arv = avrU[i]
            var sumOfDiffSq = 0.0
            arrForOneMeshPoint.forEach { real ->
                sumOfDiffSq += Math.pow(arv - real.mod, 2.0)
            }
            skoU[i] = Math.sqrt( sumOfDiffSq / size )
        }
        track.histI.forEachIndexed { i, arrForOneMeshPoint ->
            var RMS = 0.0
            val size = arrForOneMeshPoint.size
            if (rmsWindow >= size)
                throw Exception("Размер окна расчета RMS меньше количества проведенных расчетов мгновенных схем")
            for ( i in 0..rmsWindow){
                RMS += Math.pow( arrForOneMeshPoint[i].mod, 2.0 )
            }
            var maxRMS = RMS
            for (i in rmsWindow+1 .. size-1){
                RMS = RMS - Math.pow( arrForOneMeshPoint[i - rmsWindow].mod, 2.0 ) + Math.pow( arrForOneMeshPoint[i].mod, 2.0 )
                if (RMS > maxRMS) {
                    maxRMS = RMS
                }
            }
            rmsI[i] = Math.sqrt( maxRMS )
        }
    }
}

/**
 * Класс собирающий статистику (проводящий пост обрбаботку данных и потом хранит их в своих свойствах
 */
class Data (tracks: Array<Track> ) {
    internal val data = Array(tracks.size){ i -> TrackStat(tracks[i], 5) }

    fun print(track: Int){
        println("track${track}.maxU: ${Arrays.deepToString( data[track].maxU)} " )
        println("track${track}.avrU: ${Arrays.deepToString( data[track].avrU)} " )
        println("track${track}.posArvU: ${Arrays.deepToString( data[track].posAvrU)} " )
        println("track${track}.mimU: ${Arrays.deepToString( data[track].minU)} " )
        println("track${track}.skoU: ${Arrays.deepToString( data[track].skoU)} " )
        println("track${track}.rmsI: ${Arrays.deepToString( data[track].rmsI)} " )
    }
}




