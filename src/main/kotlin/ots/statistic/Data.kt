package ots.statistic

import ots.calc.Track
import java.io.File
import kotlin.math.absoluteValue
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
 * @property posAvrU массив положительных средних арифметических значений напряжений в узлах сетки, полезно для оценки потенциальной опасности электрокоррозии рельсов и железобетонных контрукций
 * @property exContactU доля времени (в %) когда было приевышено максимальное напряжения прикосновения к рельсу
 * @property maxSU массив максимальных значений напряжений статистического корридора в узлах сетки, значение напряжения попавшие в диапаон (avrU, avrU+1.5*skoU]
 * @property minSU массив минимальны значений напряжений статистического корридора в узлах сетки, значение напряжения попавшие в диапаон (avrU, avrU-1.5*skoU]
 * @property maxI массив максимальных значений тока в каждом узле сетки
 * @property rmsI массив действующих значение тока за rmsWindow шагов расчета в каждом узле сетки
 */

@Suppress("NAME_SHADOWING")
class TrackStat(track: Track, rmsWindow: Int) {
    internal val size = track.mesh.size
    internal val maxU: Array<Double> = Array(size) {0.0}
    internal val minU: Array<Double> = Array(size) {0.0}
    internal val avrU: Array<Double> = Array(size) {0.0}
    internal val skoU: Array<Double> = Array(size) {0.0}
    internal val posAvrU: Array<Double> = Array(size) {0.0}
    internal val exContactU: Array<Double> = Array(size) {0.0}
    internal val maxSU: Array<Double> = Array(size) {0.0}
    internal val minSU: Array<Double> = Array(size) {0.0}
    internal val rmsI: Array<Double> = Array(size) {0.0}
    internal val maxI: Array<Double> = Array(size) {0.0}
    /**
     * Допустимая величина потенцила рельсов,
     * равно допустимое напряжение прикосновения на теле человека (75 В) деленоек на коэффициент прикосновения (примерно 0.5)
     */
    private val railMaxAllowedU: Double = 75.0/0.5
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
            var moreThenAllowedU = 0
            arrForOneMeshPoint.forEach{ real ->
                val mod = real.mod
                if (mod < min)
                    min = mod
                if (mod > max)
                    max = mod
                sumOfUArray += mod
                if (mod > 0)
                    sumOfPosUArray += mod
                if ( mod.absoluteValue > railMaxAllowedU ){
                    moreThenAllowedU += 1
                }
            }

            maxU[i] = max
            minU[i] = min
            avrU[i] = sumOfUArray / size
            posAvrU[i] = sumOfPosUArray / size
            exContactU[i] = 100.0 * moreThenAllowedU / size
            val arv = avrU[i]
            var sumOfDiffSq = 0.0
            arrForOneMeshPoint.forEach { real ->
                sumOfDiffSq += (arv - real.mod).pow(2.0)
            }
            skoU[i] = sqrt( sumOfDiffSq / size )
            maxSU[i] = avrU[i] + 1.5 * skoU[i]
            minSU[i] = avrU[i] - 1.5 * skoU[i]
        }
        track.histI.forEachIndexed { i, arrForOneMeshPoint ->

            maxI[i] = arrForOneMeshPoint[0].mod.absoluteValue
            var sumOfSquare = 0.0
            arrForOneMeshPoint.forEach{ real ->
                val mod = real.mod.absoluteValue
                if (mod > maxI[i]){
                    maxI[i] = mod
                }
                sumOfSquare += mod.pow(2.0)
            }
            rmsI[i] = sqrt( sumOfSquare/arrForOneMeshPoint.size)
        /*
            This RMS windows algo gives rmsI value more then maxI
            maxI[i] = arrForOneMeshPoint[0].mod.absoluteValue
            var rms = 0.0
            val size = arrForOneMeshPoint.size
            if (rmsWindow >= size)
                throw Exception("Размер окна расчета RMS меньше количества проведенных расчетов мгновенных схем")
            for ( ii in 0..rmsWindow){
                val mod = arrForOneMeshPoint[ii].mod
                rms += mod.pow(2.0)
                if (mod.absoluteValue > maxI[i]){
                    maxI[i] = mod.absoluteValue
                }
            }
            var maxRMS = rms
            for ( ii in rmsWindow+1 until size){
                val mod = arrForOneMeshPoint[ii].mod
                rms = rms - arrForOneMeshPoint[ii - rmsWindow].mod.pow(2.0) + mod.pow(2.0)
                if (rms > maxRMS) {
                    maxRMS = rms
                }
                if (mod.absoluteValue > maxI[i]){
                    maxI[i] = mod.absoluteValue
                }
            }
            rmsI[i] = sqrt( maxRMS/rmsWindow )
        */
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
        println("track${track}.maxI: ${data[track].maxI.contentDeepToString()} " )
        println("track${track}.rmsI: ${data[track].rmsI.contentDeepToString()} " )
    }

    /**
     * Функция сохраняет статистику по заданному пути в указанный файл в фомате csv.
     * Каждая строка в файле содержит данные по одному узлу сетки
     */
      fun write2csv( filePath: String, trackIdx: Int ){
         File(filePath).printWriter().use { out->
             out.println("meshNode,maxU,minU,arvU,skoU,posAvrU,exContactU,maxSU,minSU,maxI,rmsI")
             val tr = data[trackIdx]
             for (i in 0 until tr.size){
                 out.println("${i},${tr.maxU[i]},${tr.minU[i]},${tr.avrU[i]},${tr.skoU[i]},${tr.posAvrU[i]},${tr.exContactU[i]},${tr.maxSU[i]},${tr.minSU[i]},${tr.maxI[i]},${tr.rmsI[i]}")
             }
         }
         println("StatData for track${trackIdx} was saved to file $filePath")
      }
}




