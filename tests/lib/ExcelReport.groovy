import java.util.zip.ZipFile
import javax.xml.parsers.DocumentBuilderFactory

/**
 * The cells of a report (an .xlsx file, a zip of XML files) as text, row by row: what the Excel file shows,
 * without its colours. For the pipeline test (tests/main.nf.test), with Java's own zip and XML readers.
 */
class ExcelReport {

    static List<List<String>> rows(String path) {
        def zip = new ZipFile(path)
        try {
            def shared = []
            def strings = zip.getEntry("xl/sharedStrings.xml")
            if (strings != null) {
                def items = parse(zip, strings).getElementsByTagName("si")
                for (int i = 0; i < items.length; i++) {
                    shared << items.item(i).textContent
                }
            }
            def out = []
            def rows = parse(zip, zip.getEntry("xl/worksheets/sheet1.xml")).getElementsByTagName("row")
            for (int r = 0; r < rows.length; r++) {
                def cells = rows.item(r).getElementsByTagName("c")
                def line = []
                for (int c = 0; c < cells.length; c++) {
                    def cell = cells.item(c)
                    int column = column(cell.getAttribute("r"))
                    while (line.size() < column) {
                        line << ""
                    }
                    def value = cell.getElementsByTagName("v")
                    def type = cell.getAttribute("t")
                    if (type == "s") {
                        line << shared[value.item(0).textContent as int]
                    } else if (type == "inlineStr") {           // an empty cell has no text in it
                        def text = cell.getElementsByTagName("is")
                        line << (text.length ? text.item(0).textContent : "")
                    } else {
                        line << (value.length ? value.item(0).textContent : "")
                    }
                }
                out << line
            }
            int width = out.collect { it.size() }.max() ?: 0
            out.each { line -> while (line.size() < width) { line << "" } }
            return out
        } finally {
            zip.close()
        }
    }

    private static parse(ZipFile zip, entry) {
        DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(zip.getInputStream(entry))
    }

    /** "C12" -> 2: the cell's column, from 0 */
    private static int column(String reference) {
        int n = 0
        for (char ch : reference.toCharArray()) {
            if (!Character.isLetter(ch)) {
                break
            }
            n = n * 26 + ((int) ch - (int) ('A' as char) + 1)
        }
        return n - 1
    }
}
