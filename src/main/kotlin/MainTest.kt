
class Mesh(
    val edge:Int,
    val dist:Double,

){
    init {
        println("create mesh obj")
    }
}

class Track(
    val name: String,
    val mesh: Mesh
){
    init {
        println("Create Track $name, mash.edge ${mesh.edge}, ${mesh.dist}")
    }
}

class Calc(
    val track: Array<Track>
)
{
    init {
        println("new trach")
        for (tr in track){
            println("${tr.name} -> mesh ${tr.mesh.edge} ${tr.mesh.dist}")
        }
    }
}

fun main(){
    println("main")
    val mesh = Mesh(1, 12.2)
    if (1>0){
        val track1 = Track("new221", mesh)
        val track2 = Track("new222", mesh)
        val calc = Calc (arrayOf(track1, track2))
    }
    val track1 = Track("new1", mesh)
    val track2 = Track("new2", mesh)
    val calc = Calc (arrayOf(track1, track2))

}