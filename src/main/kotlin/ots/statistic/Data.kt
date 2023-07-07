package ots.statistic

import ots.calc.Track
import java.io.File
import java.util.*

/**
 * Класс со статистикой по одному пути, данные распределены по узлам сетки пути
 * Исторические данные на основании которых собирается статистика храняться в класе Track
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
 * Класс собирающий статистику
 * Сбор статистики происходит в момент создания обхекта, в контрукторах TrackStat
 * @param tracks Массив путей по которым нужно собрать статистику
 */
class Data (tracks: Array<Track> ) {
    internal val data = Array(tracks.size){ i -> TrackStat(tracks[i], 5) }

    /**
     * Функция распечатывает статистку по указанному пути
     * @param track индекс пути в массиве путей
     */
    fun print(track: Int){
        println("track${track}.maxU: ${Arrays.deepToString( data[track].maxU)} " )
        println("track${track}.avrU: ${Arrays.deepToString( data[track].avrU)} " )
        println("track${track}.posArvU: ${Arrays.deepToString( data[track].posAvrU)} " )
        println("track${track}.mimU: ${Arrays.deepToString( data[track].minU)} " )
        println("track${track}.skoU: ${Arrays.deepToString( data[track].skoU)} " )
        println("track${track}.rmsI: ${Arrays.deepToString( data[track].rmsI)} " )
    }

    /**
     * Функция сохраняет статистику по заданному пути в указанный файл в фомате csv.
     * Каждая строка в файле содержит данные по одному узлу сетки
     */
      fun write2csv( filePath: String, trackIdx: Int ){
         File(filePath).printWriter().use { out->
             out.println("meshNode,maxU,minU,arvU,skoU,posAvrU,rmsI")
             val tr = data.get(trackIdx)
             for (i in 0..tr.size-1){
                 out.println("${i},${tr.maxU[i]},${tr.minU[i]},${tr.avrU[i]},${tr.skoU[i]},${tr.posAvrU[i]},${tr.rmsI[i]}")
             }
         }
     }
}




