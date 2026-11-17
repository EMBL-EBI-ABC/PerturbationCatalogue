import groovy.json.JsonOutput
import groovy.json.JsonSlurper

/** Local policy calls; computational commands remain Nextflow tasks. */
class MavedbCache {
    static Process acquire(Map engine, Object root) {
        def process = new ProcessBuilder(engine.python.toString(), engine.runtime.toString(), 'lock', '--root', root.toString())
            .redirectError(ProcessBuilder.Redirect.INHERIT).start()
        if (process.inputStream.newReader('UTF-8').readLine() != 'LOCKED') {
            process.waitFor()
            throw new IllegalStateException('Cannot acquire MaveDB output writer lock')
        }
        return process
    }

    static Map call(Map engine, List arguments, Map payload = [:]) {
        def command = [engine.python.toString(), engine.runtime.toString()] + arguments.collect { it.toString() }
        def process = new ProcessBuilder(command).redirectError(ProcessBuilder.Redirect.INHERIT).start()
        process.outputStream.withWriter('UTF-8') { writer -> writer << JsonOutput.toJson(payload) }
        def response = process.inputStream.getText('UTF-8')
        if (process.waitFor() != 0) {
            throw new IllegalStateException('MaveDB result validation failed; see the diagnostic above')
        }
        return (Map) new JsonSlurper().parseText(response)
    }

    static Map resolve(Map job) {
        def arguments = ['resolve', '--request', job.engine.request, '--stage', job.stage]
        if (job.index != null) arguments += ['--index', job.index]
        return call(job.engine, arguments, job.inputs) + [engine: job.engine, inputs: job.inputs, index: job.index]
    }

    static Map request(Map job) {
        return (Map) new JsonSlurper().parseText(new File(job.engine.request.toString()).getText('UTF-8'))
    }

    static Map specification(Map job) {
        return (Map) new JsonSlurper().parseText(new File(job.spec.toString()).getText('UTF-8'))
    }

    static Map result(Object path) {
        return (Map) new JsonSlurper().parseText(path.text)
    }

    static String arguments(List values) {
        return values.collect { quote(it) }.join(' ')
    }

    static String quote(Object value) {
        return "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }
}
