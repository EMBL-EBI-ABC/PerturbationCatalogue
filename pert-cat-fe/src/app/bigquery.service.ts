import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';

@Injectable({
  providedIn: 'root',
})
export class BigQueryService {
  private readonly bigQueryUrl = '...';

  constructor(private http: HttpClient) { }

  getData(): Observable<any> {
    const query = {
      query: `SELECT * FROM \`mavedb.metadata\` LIMIT 10`,
      useLegacySql: false
    };

    const headers = new HttpHeaders().set('Content-Type', 'application/json');

    return this.http.post(this.bigQueryUrl, query, { headers });
  }
}
