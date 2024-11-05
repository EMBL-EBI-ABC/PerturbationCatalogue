import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';


@Injectable({
  providedIn: 'root',
})
export class ElasticService {
  private readonly elasticUrl = 'https://USERNAME:KEY@SUBDOMAIN.elastic-cloud.com/mavedb/_search?filter_path=hits.hits._source';

  constructor(private http: HttpClient) { }

  getData(): Observable<any> {

    // Not working for now.
    // const query = '{"size": 10, "query": {"match_all": {}}}'
    // const headers = new HttpHeaders().set('Content-Type', 'application/json');
    // const response = this.http.post(this.elasticUrl, query, { headers })

    // Temporary substitution with first 10 entries from GitHub.
    const response = this.http.get("https://raw.githubusercontent.com/tskir/playground/refs/heads/main/test10.json")

    // Return the fields we need.
    return response.pipe(
      map((response: any) => response.hits.hits.map((hit: { _source: any; }) => hit._source))
    );
  }
}
