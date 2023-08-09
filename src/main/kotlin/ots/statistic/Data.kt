package ots.statistic

import ots.calc.Track
import java.io.File
import kotlin.math.pow
import kotlin.math.sqrt

/**
 * Класс со статистикой по одному пути, данные распределены по узлам сетки пути
 * Исторические данные на основании которых собирается статистика храняться в класе Track
 * @param track массив путей
 * @param rmsWindow количество шагов расчета за который расчитывать RMS
 * @property maxU массив максимальных значений напряжений в узлах сетки
 * @property minU массив минимальны значений напряжений в узлах сетки
 * @property avrU массив средних арифметических значений напряжений в узлах сетки
 * @property skoU массив среднеквадратичных отклонений от среднего значения напряжений у узлах сетки
 * @property posAvrU массив положительных средних арифметических значений напряжений в узлах сетки
 * @property rmsI действующее значение тока за rms_steps шагов расчета в каждом узле сетки
 */

@Suppress("NAME_SHADOWING")
class TrackStat(track: Track, rmsWindow: Int) {
    internal val size = track.mesh.size
    internal val maxU: Array<Double> = Array(size) {0.0}
    internal val minU: Array<Double> = Array(size) {0.0}
    internal val avrU: Array<Double> = Array(size) {0.0}
    internal val skoU: Array<Double> = Array(size) {0.0}
    internal val posAvrU: Array<Double> = Array(size) {0.0}
    internal val rmsI: Array<Double> = Array(size) {0.0}

    /**
     * В методе инициализации проходит вся магия приготовления/расчета статистики
     */
    init {
        track.histU.forEachIndexed { i, arrForOneMeshPoint ->
            var sumOfUArray = 0.0
            var sumOfPosUArray = 0.0
            var min = arrForOneMeshPoint[0].mod
            var max = arrForOneMeshPoint[0].mod
            val size = arrForOneMeshPoint.size
            arrForOneMeshPoint.forEach{ real ->
                val mod = real.mod
                if (mod < min)
                    min = mod
                if (mod > max)
                    max = mod
                sumOfUArray += mod
                if (mod > 0)
                    sumOfPosUArray += mod

            }
            maxU[i] = max
            minU[i] = min
            avrU[i] = sumOfUArray / size
            posAvrU[i] = sumOfPosUArray / size
            val arv = avrU[i]
            var sumOfDiffSq = 0.0
            arrForOneMeshPoint.forEach { real ->
                sumOfDiffSq += (arv - real.mod).pow(2.0)
            }
            skoU[i] = sqrt( sumOfDiffSq / size )
        }
        track.histI.forEachIndexed { i, arrForOneMeshPoint ->
            var rms = 0.0
            val size = arrForOneMeshPoint.size
            if (rmsWindow >= size)
                throw Exception("Размер окна расчета RMS меньше количества проведенных расчетов мгновенных схем")
            for ( i in 0..rmsWindow){
                rms += arrForOneMeshPoint[i].mod.pow(2.0)
            }
            var maxRMS = rms
            for ( i in rmsWindow+1 until size){
                rms = rms - arrForOneMeshPoint[i - rmsWindow].mod.pow(2.0) + arrForOneMeshPoint[i].mod.pow(2.0)
                if (rms > maxRMS) {
                    maxRMS = rms
                }
            }
            rmsI[i] = sqrt( maxRMS )
        }
    }
}

/**
 * Класс собирающий статистику
 * Сбор статистики происходит в момент создания объекта, в конструкторах TrackStat
 * @param tracks Массив путей по которым нужно собрать статистику
 */
class Data (tracks: Array<Track> ) {
    private val data = Array(tracks.size){ i -> TrackStat(tracks[i], 5) }

    /**
     * Функция распечатывает статистку по указанному пути
     * @param track индекс пути в массиве путей
     */
    fun print(track: Int){
        println("track${track}.maxU: ${data[track].maxU.contentDeepToString()} " )
        println("track${track}.avrU: ${data[track].avrU.contentDeepToString()} " )
        println("track${track}.posArvU: ${data[track].posAvrU.contentDeepToString()} " )
        println("track${track}.mimU: ${data[track].minU.contentDeepToString()} " )
        println("track${track}.skoU: ${data[track].skoU.contentDeepToString()} " )
        println("track${track}.rmsI: ${data[track].rmsI.contentDeepToString()} " )
    }

    /**
     * Функция сохраняет статистику по заданному пути в указанный файл в фомате csv.
     * Каждая строка в файле содержит данные по одному узлу сетки
     */
      fun write2csv( filePath: String, trackIdx: Int ){
         File(filePath).printWriter().use { out->
             out.println("meshNode,maxU,minU,arvU,skoU,posAvrU,rmsI")
             val tr = data[trackIdx]
             for (i in 0 until tr.size){
                 out.println("${i},${tr.maxU[i]},${tr.minU[i]},${tr.avrU[i]},${tr.skoU[i]},${tr.posAvrU[i]},${tr.rmsI[i]}")
             }
         }
         println("StatData for track${trackIdx} was saved to file $filePath")
      }
}




